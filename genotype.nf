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
This workflow takes bam files as input and outputs a multi-sample vcf and related stats.
The multi-sample vcf can be used as input to BQSR workflow to save running GenomicsDBImport and GenotypeGVCFs
as the first two processes of the BQSR pipeline.
The multi-sample vcf is the vcfs that will be carried forwards for analysis.
This workflow outputs summary metrics of the produced multi-sample vcf file that can be used to compare the effect
running BQSR (e.g. by comparing Ti/Tv ratios or number of SNPs).
*/

bin_dir = "${workflow.launchDir}/bin"
params.split_haplotype_caller = false

params.overwrite = false

params.output_dir = "${workflow.launchDir}/outputs/gatk_variant_calling"
gatk_output_vcf_publishDir = [params.output_dir, "gatk_output_vcf"].join(File.separator)
gatk_vcf_stats_publishDir = [params.output_dir, "gatk_vcf_stats"].join(File.separator)
bcftools_call_publishDir = [params.output_dir, "bcftools_output_vcf"].join(File.separator)
bcftools_vcf_stats_publishDir = [params.output_dir, "bcftools_vcf_stats"].join(File.separator)

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

// Create a channel that is a tuple of the pair_id and the corresponding bam as output from the preprocessing workflow
Channel.fromFilePairs("${params.input_dir}/*.bam{,.bai}").toList().flatMap{
        list_of_bams -> 
        
        def group_list = [];
        // For each of the bam sets, pull out the current group name and put it into the group_list
        list_of_bams.eachWithIndex{ bam_set, i -> group_list << bam_set[0];}
        
        // Then work backwards on a per character basis checking to see if the character is found in common
        // for all samples. If it is then this will be a character to discard. Keep going until
        // we find a character that is not in common across all names and make a note of the index position.
        // We will then return the group ids with this number of characters discarded
        def cut_index = 0;
        def exit = false;
        for (i = 1; !exit; i++) {
            // For each character starting from the last
            start_char = group_list[0].charAt(group_list[0].length() - i)
            
            group_list.eachWithIndex{
                group_name, j -> 
                    if (group_name.charAt(group_name.length() - i) != start_char){
                        // Then j is the index to start cutting from
                        cut_index = -1 * i;
                        
                        exit = true;
                    }
                }
        }
        // Now that we have the cut index we can cut each of the group names
        // At this point we want to create a new list and return it
        def output_list = [];
        list_of_bams.eachWithIndex{bam_set, i -> output_list << [bam_set[0][0..cut_index], bam_set[1] ];}
        return output_list;
    }.set{input_bam_ch}

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


process index_dictionary_refgenome{
    container 'broadinstitute/gatk:4.2.0.0'

    input:
    path ref_genome from params.ref

    output:
    tuple path(ref_genome), path("*.dict"), path("*.fai") into index_dictionary_refgenome_out_ch

    script:
    """
    gatk CreateSequenceDictionary -R ${ref_genome}
    samtools faidx ${ref_genome}
    """
}

// Control the forking of the genotyping by forking the pipelines
if (params.mode == "bcftools"){
    index_dictionary_refgenome_out_ch.into{bcftools_mpileup_ref_genome_ch; bcftools_vcfstats_ref_genome_ch}
    input_bam_ch.set{bcftools_mpileup_ch}
}else if (params.mode == "gatk"){
    index_dictionary_refgenome_out_ch.into{gatk_haplotype_caller_ref_genome_ch; GenotypeGVCFs_ref_genome_ch; bcftools_vcfstats_ref_genome_ch}
    input_bam_ch.set{gatk_haplotype_caller_gvcf_ch}
}else{
    // Default = "both"
    params.mode = "both"
    index_dictionary_refgenome_out_ch.into{gatk_haplotype_caller_ref_genome_ch; bcftools_mpileup_ref_genome_ch; GenotypeGVCFs_ref_genome_ch; bcftools_vcfstats_ref_genome_ch}
    input_bam_ch.into{bcftools_mpileup_ch; gatk_haplotype_caller_gvcf_ch}
}
if (params.mode == "both" || params.mode == "bcftools"){
    // START OF BCFTOOLS MPILEUP AND CALLING
    // The bcftools mode will also follow a scatter-gather approach
    // For each scaffold of the reference run bcftools mpileup
    // Then combine with bcftools concat
    // finally run bcftools call
    process bcftools_mpileup{
        container 'halllab/bcftools:v1.9'

        input:
        each scaffold from Channel.fromList(scaffold_list)
        tuple path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from bcftools_mpileup_ref_genome_ch
        path(bam_bai) from bcftools_mpileup_ch.collect{it[1]}

        output:
        tuple val(scaffold), path("${scaffold}.vcf.gz"), path("${scaffold}.vcf.gz.csi") into bcftools_call_ch

        script:
        """
        bcftools mpileup -f ${ref_genome} -d 100 -Oz -r $scaffold -o ${scaffold}.vcf.gz ${bam_bai.findResults { i -> if(i.getName().endsWith(".bam")){"$i "}else{null}}.join()}
        bcftools index ${scaffold}.vcf.gz
        """
    }

    process bcftools_call{
        container 'halllab/bcftools:v1.9'
        
        input:
        tuple val(scaffold), path(scaffold_vcf), path(scaffold_vcf_csi) from bcftools_call_ch

        output:
        path("${scaffold}.call.vcf.gz") into bcftools_concat_ch
        path("${scaffold}.call.vcf.gz.csi") into bcftools_concat_csi_ch

        script:
        """
        bcftools call -Oz -o ${scaffold}.call.vcf.gz -r ${scaffold} -m -v $scaffold_vcf
        bcftools index ${scaffold}.call.vcf.gz
        """

    }

    process bcftools_concat{
        container 'halllab/bcftools:v1.9'
        publishDir bcftools_call_publishDir, mode: "copy"

        input:
        path(per_scaffold_bcfs) from bcftools_concat_ch.collect()
        path(per_scaffold_bcfs_csi) from bcftools_concat_csi_ch.collect()

        output:
        tuple path("bcftools.call.concat.vcf.gz"), path("bcftools.call.concat.vcf.gz.tbi") into bcftools_out_ch

        script:
        """
        bcftools concat -Oz -o bcftools.call.concat.vcf.gz ${per_scaffold_bcfs.collect{ i -> "$i "}.join()}
        tabix -p vcf bcftools.call.concat.vcf.gz
        """
    }
}

if (params.mode == "both" || params.mode == "gatk"){
    // START OF GATK HAPLOTYPE CALLING
    // Call variants on a per scaffold (chromosone) basis and then gather back to a single vcf (scatter-gather)
    // 1 - Call haplotypes on a per sample basis
    // 2 - run GenomicsDBImport on a per scaffold basis
    // 3 - run GenotypeGVCFs on a per scaffold basis.
    // 4 - Gather the per scaffold files using GatherVcfs
    // 5 - Evaluate the vcfs using bcftools stats and rtg vcfstats
    // NB Good reference for the scatter-gather strategy employed here: https://github.com/IARCbioinfo/gatk4-GenotypeGVCFs-nf/blob/master/gatk4-GenotypeGVCFs.nf
    // NB HaplotypCaller requires a .fai
    // https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
    if (params.split_haplotype_caller){
        process gatk_haplotype_caller_gvcf_split{
            tag {pair_id}
            container 'broadinstitute/gatk:4.2.0.0'
            cpus Math.max(params.gatk_haplotype_caller_cpus - 1 , 1)

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

    // NB Multiple cores seem to be used by each instance of the process
    // I cannot see a way to control the thread usage via the command line options
    // As such, I have upped the required number of cpus for this process to 3
    process GenotypeGVCFs{
        container 'broadinstitute/gatk:4.2.0.0'
        cpus 3
        tag "${scaffold}"

        input:
        tuple val(scaffold), file(workspace), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from genotype_GVCFs_ch.combine(GenotypeGVCFs_ref_genome_ch)

        output:
        file("GenotypeGVCFs.out.${scaffold}.vcf") into gather_vcfs_for_eval_ch
        file("GenotypeGVCFs.out.${scaffold}.vcf.idx") into gather_vcfs_idx_for_eval_ch

        script:
        """
        WORKSPACE=\$( basename ${workspace} )
        gatk GenotypeGVCFs -R ${ref_genome} -O GenotypeGVCFs.out.${scaffold}.vcf \
        --only-output-calls-starting-in-intervals -V gendb://\$WORKSPACE -L ${scaffold}
        """
    }

    process gather_vcfs_for_eval{
        tag "GatherVcfs_for_eval"
        container 'broadinstitute/gatk:4.2.0.0'
        publishDir gatk_output_vcf_publishDir, mode: 'copy', pattern: '*.{gz,gz.tbi}'

        input:
        val(scaffolds) from Channel.fromList(scaffold_list).collect()
        path(vcf) from gather_vcfs_for_eval_ch.collect()
        path(vcf_idx) from gather_vcfs_idx_for_eval_ch.collect()

        output:
        tuple path("gatk.out.vcf.gz"), path("gatk.out.vcf.gz.tbi") into gatk_out_ch

        script:
        """
        gatk GatherVcfs ${scaffolds.collect{ "--INPUT GenotypeGVCFs.out.${it}.vcf " }.join()} --OUTPUT gatk.out.vcf
        bgzip -c gatk.out.vcf > gatk.out.vcf.gz
        tabix -p vcf gatk.out.vcf.gz
        """
    }
}

if (params.mode == "bcftools"){
    bcftools_out_ch.into{bcftools_stats_ch; rtg_vcfstats_ch}
}else if (params.mode == "gatk"){
    gatk_out_ch.into{bcftools_stats_ch; rtg_vcfstats_ch}
}else if (params.mode == "both"){
    // Default
    bcftools_out_ch.concat(gatk_out_ch).into{bcftools_stats_ch; rtg_vcfstats_ch}
}

// VCF STATISTICS
// useful reference for rtg: https://genomics.sschmeier.com/ngs-variantcalling/index.html
process bcftools_vcfstats{
    tag "bcftools_stats"
    container "halllab/bcftools:v1.9"
    publishDir gatk_vcf_stats_publishDir, mode: "copy", pattern: "*gatk*"
    publishDir bcftools_vcf_stats_publishDir, mode: "copy", pattern: "*bcftools.call*"

    input:
    tuple path(vcfgz), path(vcfgz_tbi), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from bcftools_stats_ch.combine(bcftools_vcfstats_ref_genome_ch)

    output:
    path("*.bcftools.stats.txt") into bcftools_stats_out_ch

    script:
    stats_out_file = vcfgz.getName().replaceAll(".vcf.gz", ".bcftools.stats.txt")
    """
    bcftools stats -F $ref_genome -s - $vcfgz > $stats_out_file
    """
}

// This produces a really handy per sample output but ideally we'd also like a summary metric (i.e. across all samples)
// We have written and use summarise_rtg_vcfstats.py to do this.
// It outputs a file rtg.stats.summary*.txt (with optional iteration annotation) that contains the total
// SNPs and weighted averges (SNPs per sample) of a selection of metrics including Ti/Tv ratio and insertion/deletion ratio
process rtg_vcfstats_per_sample{
    tag "rtg_vcfstats"
    container 'realtimegenomics/rtg-tools:latest'
    publishDir gatk_vcf_stats_publishDir, mode: "copy", pattern: "*gatk*"
    publishDir bcftools_vcf_stats_publishDir, mode: "copy", pattern: "*bcftools.call*"

    input:
    tuple path(vcfgz), path(vcfgz_tbi) from rtg_vcfstats_ch

    output:
    path("*.rtg.stats.per_sample.txt") into rtg_vcfstats_out_ch

    script:
    stats_out_file = vcfgz.getName().replaceAll(".vcf.gz", ".rtg.stats.per_sample.txt")
    """
    rtg vcfstats $vcfgz > $stats_out_file
    """
}

process rtg_vcfstats_summary{
    tag "rtg_vcfstats_summary"
    container 'broadinstitute/gatk:4.2.0.0'
    publishDir gatk_vcf_stats_publishDir, mode: "copy", pattern: "*gatk*"
    publishDir bcftools_vcf_stats_publishDir, mode: "copy", pattern: "*bcftools.call*"

    input:
    path(rtg_vcfstats_output) from rtg_vcfstats_out_ch

    output:
    path("*rtg.stats.summary.txt") into rtg_vcfstats_summary_out_ch

    script:
    stats_out_file = rtg_vcfstats_output.getName().replaceAll("per_sample", "summary")
    """
    python3 ${bin_dir}/summarise_rtg_vcfstats.py $rtg_vcfstats_output $stats_out_file
    """
}
