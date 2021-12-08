#!/usr/bin/env nextflow

/*
    Pipeline to assess relatedness of non-model organism samples through
    shallow whole genome sequencing.

    Copyright (C) 2021  Benjamin C. C. Hume

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Author: Benjamin C. C. Hume
    Contact: didillysquat@gmail.com
*/

/*
This workflow is used to perform a single round of BQSR. BQSR corrects aligned sequencing data (i.e. bam files)
according to a number of covariates. To do this it needs to have a set of known variants as input.
Because there are not a reference set of variants 
available, a set of 'high confidence' variants are generated by hard filtering and used as the 'known' variants.
The BQSR process is in two steps, the first generating a recalibration table and the second applying this table
to the existing bams to create a recalibrated set of bams.

As input, this workflow requires the set of .bam files that BQSR is to be performed on. The input directory should contain
paired .bam and .bam.bai files. As an optional input, a multi-sample vcf can be provided upon which the hardfiltering
will be performed to generate the 'known variants'. If this vcf is not provided, the first three stages of the workflow
will produce a multi-sample vcf on a per scaffold basis. This part of the process is computationally intensive. As
such, where possible the vcf should be provided. If this workflow is being run after the 'gatk call variants' workflow
the vcf output from this process may be provided. The input should point directly to the .vcf file and the .vcf.idx
should be present in the same directory.

This workflow outputs a set of .bam files that have had BSR applied to them. It also outputs a set of metrics detailing
the differences before and after the BQSR.
*/

bin_dir = "${workflow.launchDir}/bin"
params.overwrite = false
params.split_haplotype_caller = false

params.output_dir = "${workflow.launchDir}/outputs/bqsr"
gatk_bqsr_vcf_publishDir = [params.output_dir, "gatk_bqsr_output_vcf"].join(File.separator)
gatk_bqsr_bam_publishDir = [params.output_dir, "gatk_bqsr_output_bams"].join(File.separator)
gatk_bqsr_vcf_stats_publishDir = [params.output_dir, "gatk_bqsr_vcf_stats"].join(File.separator)
gatk_bqsr_analyze_covariates_publishDir = [params.output_dir, "gatk_bqsr_analyze_covariates_metrics"].join(File.separator)

if (!params.haplotypecaller_max_mem){
    params.haplotypecaller_max_mem = 8
}
if (!params.gatk_haplotype_caller_cpus){
    params.gatk_haplotype_caller_cpus = 4
}

// We will enforce a check to make sure that the output directory doesn't already exist as we don't want to
// accidentally overwrite the files. We will only do this check if overwrite is false.
if(!params.overwrite){
    // Check to see if each of the output dirs already exist
    def path_to_check = new File(params.output_dir);
    if (path_to_check.exists()) {
        throw new Exception("""The output directory ${path_to_check} already exists.\n
        If you want to overwrite the contents of this directory, please provide the --overwrite flag.\n""")
    }
}

// This scaffold list is created from the reference fasta and is used for the scatter-gather approach
// used in vcf calling.
scaffold_list = {  
    def scaffolds = []
    new File(params.ref).eachLine {
        line -> 
        if (line.startsWith(">")){scaffolds << line.split()[0][1..-1];}
        }
    return scaffolds
}()

// Call variants on a per scaffold (chromosone) basis.
// Follow a scatter-gather strategy 
// (first three steps can be skipped if provided with a vcf to hard filter from)
// (if a vcf is provided at input, an additional step 4a splits the vcf per scaffold for filtering)
// 1 - Call haplotypes on a per sample basis
// 2 - run GenomicsDBImport on a per scaffold basis
// 3 - run GenotypeGVCFs on a per scaffold basis.
// 4a - (if working from provided vcf) split vcf by scaffold
// 4 - Hard filter the vcfs on a per scaffold basis using VariantFiltration and remove the filtered variants using SelectVariants
// 5 - Gather the per scaffold files using GatherVcfs GatherVcfs to collect the per chromosome files into a single vcf
// 6 - Make the BQSR recallibration table on a persample basis with BaseRecalibrator using the vcf from 5 as known variants
// 7 - Apply the recalibration table from 6 on a per sample basis with  ApplyBQSR
// 8 - Make a second recalibration table, using the new bam files from 7 as input and the output from 5 as known variants using BaseRecalibrator
// 9 - Produce metrics for the change between the bam files using before and after recalibration using the two tables as input using AnalyzeCovariates
// Good reference for the nextflow coding: https://github.com/IARCbioinfo/gatk4-GenotypeGVCFs-nf/blob/master/gatk4-GenotypeGVCFs.nf

process index_dictionary_refgenome{
    container "broadinstitute/gatk:4.2.0.0"
    
    input:
    path ref_genome from params.ref

    output:
    tuple path(ref_genome), path("*.dict"), path("*.fai") into gatk_haplotype_caller_ref_genome_ch,GenotypeGVCFs_ref_genome_ch,make_bqsr_tables_ref_genome_ch,apply_bqsr_tables_ref_genome_ch,make_bqsr_tables_ref_genome_round_2_ch,split_vcf_by_scaffold_ref_genome_ch,bcftools_vcfstats_ref_genome_ch

    script:
    """
    gatk CreateSequenceDictionary -R ${ref_genome}
    samtools faidx ${ref_genome}
    """
    }

params.input_vcf = false

// We will generate the group name by removing all parts of the file name that are found in common between all samples.
bam_ch = Channel.fromFilePairs("${params.input_dir}/*.bam{,.bai}").toList().flatMap{
    list_of_bams -> 
    // println("This is the list_of_bams: ${list_of_bams}")
    def group_list = [];
    // For each of the bam sets, pull out the current group name and put it into the group_list
    list_of_bams.eachWithIndex{ bam_set, i -> group_list << bam_set[0];}
    // println("This is the group_list: ${group_list}" )
    // Then work backwards on a per character basis checking to see if the character is found in common
    // for all samples. If it is then this will be a character to discard. Keep going until
    // we find a character that is not in common across all names and make a note of the index position.
    // We will then return the group ids with this number of characters discarded
    def cut_index = 0;
    def exit = false;
    for (i = 1; !exit; i++) {
        // For each character starting from the last
        start_char = group_list[0].charAt(group_list[0].length() - i)
        // println("Start char is currently: ${start_char}")
        group_list.eachWithIndex{
            group_name, j -> 
                if (group_name.charAt(group_name.length() - i) != start_char){
                    // Then j is the index to start cutting from
                    cut_index = -1 * i;
                    // println("Look we got here and the cut index is: ${cut_index}")
                    exit = true;
                }
            }
    }
    // Now that we have the cut index we can cut each of the group names
    // At this point we want to create a new list and return it
    def output_list = [];
    list_of_bams.eachWithIndex{bam_set, i -> output_list << [bam_set[0][0..cut_index], bam_set[1] ];}
    return output_list;
}

if (!params.input_vcf){
    // Then we are working from bams only and need to run all processes.
    bam_ch.into{gatk_haplotype_caller_gvcf_ch; make_bqsr_tables_bam_ch}
    
    
    // // NB HaplotypCaller requires a .fai
    // // https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
    // TODO try to identify a way to better parallelise this process as it is very innefficient like this
    // Work on a per scaffold basis to start with.
    if (params.split_haplotype_caller){
        process gatk_haplotype_caller_gvcf_split{
            tag {pair_id}
            container 'broadinstitute/gatk:4.2.0.0'
            cpus params.gatk_haplotype_caller_cpus

            input:
            each scaffold from Channel.fromList(scaffold_list)
            tuple val(pair_id), path(merged), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from gatk_haplotype_caller_gvcf_ch.combine(gatk_haplotype_caller_ref_genome_ch)

            output:
            tuple val(pair_id), val(scaffold), path("${pair_id}.${scaffold}.merged.g.vcf.gz") into genomics_db_import_ch
            tuple val(pair_id), val(scaffold), path("${pair_id}.${scaffold}.merged.g.vcf.gz.tbi") into genomics_db_import_tbi_ch

            script:
            """
            gatk --java-options "-Xmx${params.haplotypecaller_max_mem}g" HaplotypeCaller --native-pair-hmm-threads ${task.cpus} -R $ref_genome -I ${merged[0]} -O ${pair_id}.${scaffold}.merged.g.vcf.gz -ERC GVCF -L $scaffold
            """
        }

        process genomics_db_import_split{
            tag "${scaffold}"
            cpus 1
            container 'broadinstitute/gatk:4.2.0.0'

            input:
            tuple val(pair_id_list), val(scaffold), path(gvcf_list) from genomics_db_import_ch.groupTuple(by: 1)
            tuple val(pair_id), val(scaffold), path(gvcf_tbi) from genomics_db_import_tbi_ch.groupTuple(by: 1)
            
            output:
            tuple val(scaffold), file("genomicsdbi.out.${scaffold}") into genotype_GVCFs_ch
            
            script:
            """
            gatk  --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenomicsDBImport ${gvcf_list.collect { "-V $it " }.join()} --genomicsdb-workspace-path genomicsdbi.out.${scaffold} -L $scaffold --batch-size 50
            """
        }
    }else{
        process gatk_haplotype_caller_gvcf_no_split{
            tag {pair_id}
            container 'broadinstitute/gatk:4.2.0.0'
            cpus params.gatk_haplotype_caller_cpus

            input:
            tuple val(pair_id), path(merged), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from gatk_haplotype_caller_gvcf_ch.combine(gatk_haplotype_caller_ref_genome_ch)

            output:
            tuple val(pair_id), path("${pair_id}.merged.g.vcf.gz") into genomics_db_import_ch
            tuple val(pair_id), path("${pair_id}.merged.g.vcf.gz.tbi") into genomics_db_import_tbi_ch

            script:
            """
            gatk --java-options "-Xmx${params.haplotypecaller_max_mem}g" --native-pair-hmm-threads ${task.cpus} HaplotypeCaller -R $ref_genome -I ${merged[0]} -O ${pair_id}.merged.g.vcf.gz -ERC GVCF
            """
        }

        process genomics_db_import_no_split{
            tag "${scaffold}"
            cpus 1
            container 'broadinstitute/gatk:4.2.0.0'

            input:
            each scaffold from Channel.fromList(scaffold_list)
            path(gvcf) from genomics_db_import_ch.collect{it[1]}
            path(gvcf_tbi) from genomics_db_import_tbi_ch.collect{it[1]}
            
            output:
            tuple val(scaffold), path("genomicsdbi.out.${scaffold}") into genotype_GVCFs_ch
            
            script:
            """
            gatk  --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenomicsDBImport ${gvcf.collect { "-V $it " }.join()} --genomicsdb-workspace-path genomicsdbi.out.${scaffold} -L $scaffold --batch-size 50
            """
        }
    }

    process GenotypeGVCFs{
        container 'broadinstitute/gatk:4.2.0.0'
        cpus 5
        tag "${scaffold}"

        input:
        tuple val(scaffold), file(workspace), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from genotype_GVCFs_ch.combine(GenotypeGVCFs_ref_genome_ch)

        output:
        tuple val(scaffold), file("GenotypeGVCFs.out.${scaffold}.vcf"), file("GenotypeGVCFs.out.${scaffold}.vcf.idx") into hard_filter_ch,gather_vcf_output

        script:
        """
        WORKSPACE=\$( basename ${workspace} )
        gatk GenotypeGVCFs -R ${ref_genome} -O GenotypeGVCFs.out.${scaffold}.vcf \
        --only-output-calls-starting-in-intervals -V gendb://\$WORKSPACE -L ${scaffold}
        """
    }
    
    // NB GatherVcfs must be provided with the vcf files in order of the scaffolds. We do this using the scaffold_list.
    process gather_vcfs_for_eval{
    	tag "GatherVcfs_for_eval"
        container 'broadinstitute/gatk:4.2.0.0'
        publishDir gatk_bqsr_vcf_publishDir, mode: 'copy'

        input:
        val(scaffolds) from Channel.fromList(scaffold_list).collect()
        path(vcf_files) from gather_vcf_output.collect{[it[1], it[2]]}

    	output:
        tuple path("out.vcf"), path("out.vcf.idx") into bcftools_stats_ch,rtg_vcfstats_ch

        script:
    	"""
    	gatk GatherVcfs ${scaffolds.collect{ "--INPUT GenotypeGVCFs.out.${it}.vcf " }.join()} --OUTPUT out.vcf
    	"""
    }

    // useful reference for rtg: https://genomics.sschmeier.com/ngs-variantcalling/index.html
    process bcftools_vcfstats{
        tag "bcftools_stats"
        container "halllab/bcftools:v1.9"
        publishDir gatk_bqsr_vcf_stats_publishDir, mode: "copy"

        input:
        tuple path(vcf), path(vcf_idx), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from bcftools_stats_ch.combine(bcftools_vcfstats_ref_genome_ch)

        output:
        path("bcftools.stats.txt") into bcftools_stats_out_ch

        script:
        """
        bcftools stats -F $ref_genome -s - $vcf > bcftools.stats.txt
        """
    }

    // This produces a really handy per sample output but ideally we'd also like a summary metric
    // We have written and use summarise_rtg_vcfstats.py to do this.
    // It outputs a file rtg.stats.summary*.txt that contains the total
    // SNPs and weighted averges (SNPs per sample) of a selection of metrics including Ti/Tv ratio and insertion/deletion ratio
    process rtg_vcfstats_per_sample{
        tag "rtg_vcfstats"
        container 'realtimegenomics/rtg-tools:latest'
        publishDir gatk_bqsr_vcf_stats_publishDir, mode: "copy"

        input:
        tuple path(vcf), path(vcf_idx) from rtg_vcfstats_ch

        output:
        path("rtg.stats.per_sample.txt") into rtg_vcfstats_out_ch

        script:
        """
        rtg vcfstats $vcf > rtg.stats.per_sample.txt
        """
    }

    process rtg_vcfstats_summary{
        tag "rtg_vcfstats_summary"
        container 'broadinstitute/gatk:4.2.0.0'
        publishDir gatk_bqsr_vcf_stats_publishDir, mode: "copy"

        input:
        path(rtg_vcfstats_output) from rtg_vcfstats_out_ch

        output:
        path("rtg.stats.summary*.txt") into rtg_vcfstats_summary_out_ch

        script:
        """
        python3 ${bin_dir}/summarise_rtg_vcfstats.py $rtg_vcfstats_output
        """
    }
}else{
    // Then we have been provided with a vcf to hard filter from.
    // Can skip HaplotypeCaller, GenomicsDBImport and GenotypeGVCFs
    // Need to split the provided vcf by haplotype to pass into hard filtering
    // Output of the process needs to be a channel of tuples of the scaffold name, the vcf and the vcf.idx
    bam_ch.set{make_bqsr_tables_bam_ch}
    split_vcf_by_scaffold_vcf_ch = Channel.fromPath(params.input_vcf)
    split_vcf_by_scaffold_vcf_idx_ch = Channel.fromPath("${params.input_vcf}.idx")
    
    process split_vcf_by_scaffold{
        tag {scaffold}
        cpus 1
        memory "24GB"
        container 'broadinstitute/gatk:4.2.0.0'

        input:
        each scaffold from Channel.fromList(scaffold_list)
        path(vcf) from split_vcf_by_scaffold_vcf_ch
        tuple path(vcf_idx), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from split_vcf_by_scaffold_vcf_idx_ch.combine(split_vcf_by_scaffold_ref_genome_ch)

        output:
        tuple val(scaffold), file("GenotypeGVCFs.out.${scaffold}.vcf"), file("GenotypeGVCFs.out.${scaffold}.vcf.idx") into hard_filter_ch

        script:
        """
        gatk SelectVariants -R $ref_genome -V $vcf -L $scaffold -O "GenotypeGVCFs.out.${scaffold}.vcf"
        """

    }

}

process hard_filter{
	tag {scaffold}
    container 'broadinstitute/gatk:4.2.0.0'
    cpus 5

    input:
	tuple val(scaffold), path(vcf), path(vcfidx) from hard_filter_ch

	output:
    file("${scaffold}.filtered.vcf") into vcf_gather_vcfs_ch
    file("${scaffold}.filtered.vcf.idx") into vcf_idx_gather_vcfs_ch

    script:
	"""
	gatk VariantFiltration \
    -filter "Qual >= 100" --filter-name "Qual100" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "MQ < 35.0" --filter-name "MQ35" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "HaplotypeScore > 13.0" --filter-name "Haplo13" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRank-8" \
    -V ${vcf} -O ${scaffold}.markfiltered.vcf

    gatk SelectVariants --exclude-filtered \
      -V ${scaffold}.markfiltered.vcf \
      -O ${scaffold}.filtered.vcf
	"""
}

// // NB GatherVcfs must be provided with the vcf files in order of the scaffolds. We do this using the scaffold_list.
process gather_vcfs{
	tag "GatherVcfs"
    container 'broadinstitute/gatk:4.2.0.0'

    input:
    val(scaffolds) from Channel.fromList(scaffold_list).collect()
    path(vcf) from vcf_gather_vcfs_ch.collect()
	path(vcf_idx) from vcf_idx_gather_vcfs_ch.collect()

	output:
    tuple path("fowl.filtered.vcf"), path("fowl.filtered.vcf.idx") into make_bqsr_tables_known_variants_ch,make_bqsr_tables_known_variants_round_2_ch

    script:
	"""
	gatk GatherVcfs ${scaffolds.collect{ "--INPUT ${it}.filtered.vcf " }.join()} --OUTPUT fowl.filtered.vcf
	"""
}


// // We will want to do this on a per sample basis
// // We will want to use the fowl.vcf and the fowl.vcf.idx multiple times withhout using them up
// // It may be that we have to output these files as values rather than path/files in the gather_vcf process.
// // make_bqsr_tables_known_variants_ch.view()
// NB The BaseRecalibrator process seems to use multiple threads and uses a lot of memory
// I am struggling to get a logical response to setting the memory directive but setting the cpus to 4
// appears to be working quite well. If not set, we run out of memory
process make_bqsr_tables{
    tag {pair_id}
    container 'broadinstitute/gatk:4.2.0.0'
    cpus 4

    input:
    tuple val(pair_id), path(merged_bam), path(known_variants), path(known_variants_idx), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from make_bqsr_tables_bam_ch.combine(make_bqsr_tables_known_variants_ch).combine(make_bqsr_tables_ref_genome_ch)

    output:
    tuple val(pair_id), path(merged_bam), path("${pair_id}.table") into apply_bqsr_tables_ch
    tuple val(pair_id), path("${pair_id}.table") into compare_tables_round_1_ch

    script:
    """
    gatk BaseRecalibrator -I ${merged_bam[0]} -R $ref_genome \
    --known-sites $known_variants -O ${pair_id}.table
    """
}

process apply_bqsr_tables{
    tag {pair_id}
    container 'broadinstitute/gatk:4.2.0.0'
    publishDir gatk_bqsr_bam_publishDir, mode: "copy"

    input:
    tuple val(pair_id), path(merged_bam), path(table), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from apply_bqsr_tables_ch.combine(apply_bqsr_tables_ref_genome_ch)

    output:
    tuple val(pair_id), path("${pair_id}*.recalibrated*.bam{,.bai}") into apply_bqsr_tables_out_ch,apply_bqsr_tables_round_2_ch

    script:
    out_bam = merged_bam[0].getName().replaceAll(".bam", ".recalibrated.bam")
    out_bam_index_original = merged_bam[0].getName().replaceAll(".bam", ".recalibrated.bai")
    out_bam_index_new = merged_bam[0].getName().replaceAll(".bam", ".recalibrated.bam.bai")
    """
    gatk ApplyBQSR -R $ref_genome -I ${merged_bam[0]} \
    --bqsr-recal-file $table \
    -O ${out_bam} --create-output-bam-index
    
    mv ${out_bam_index_original} ${out_bam_index_new}
    """
}

// // See the below issue to get an idea of the general strategy of getting AnalyzeCovariates to work in 
// // GATK 4.
// // https://github.com/broadinstitute/gatk/issues/322

process make_bqsr_tables_round_2{
    tag {pair_id}
    container 'broadinstitute/gatk:4.2.0.0'
    cpus 4

    input:
    tuple val(pair_id), path(merged_bam), path(known_variants), path(known_variants_idx), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from apply_bqsr_tables_round_2_ch.combine(make_bqsr_tables_known_variants_round_2_ch).combine(make_bqsr_tables_ref_genome_round_2_ch)

    output:
    tuple val(pair_id), path("${pair_id}.after.table") into compare_tables_round_2_ch

    script:
    """
    gatk BaseRecalibrator -I ${merged_bam[0]} -R $ref_genome \
    --known-sites $known_variants -O ${pair_id}.after.table
    """
}

process AnalyzeCovariates{
    tag {pair_id}
    container 'broadinstitute/gatk:4.2.0.0'
    publishDir gatk_bqsr_analyze_covariates_publishDir, mode: "copy"

    input:
    tuple val(pair_id), path(table_before), path(table_after) from compare_tables_round_1_ch.join(compare_tables_round_2_ch)

    output:
    tuple val(pair_id), path("${pair_id}.BQSR.csv"), path("${pair_id}.BQSR.pdf") into compare_tables_round_2_out_ch
    
    script:
    """
    gatk AnalyzeCovariates \
      -before $table_before \
      -after $table_after \
      -csv ${pair_id}.BQSR.csv \
      -plots ${pair_id}.BQSR.pdf
    """
}