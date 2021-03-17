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
params.iteration = 0
params.overwrite = false
if (params.subsample){
    params.output_dir = "${workflow.launchDir}/sub_sampled_outputs/gatk_variant_calling_${params.iteration}"
    gatk_output_vcf_publishDir = [params.output_dir, "gatk_output_variants_sub_${params.subsample_depth}"].join(File.separator)
    gatk_vcf_stats_publishDir = [params.output_dir, "gatk_vcf_stats_sub_${params.subsample_depth}"].join(File.separator)
}else{
    params.output_dir = "${workflow.launchDir}/outputs/gatk_variant_calling_${params.iteration}"
    gatk_output_vcf_publishDir = [params.output_dir, "gatk_output_variants"].join(File.separator)
    gatk_vcf_stats_publishDir = [params.output_dir, "gatk_vcf_stats"].join(File.separator)
}

// We will enforce a check to make sure that the output directory doesn't already exist as we don't want to
// accidentally overwrite the files. We will only do this check if overwrite is false.
if(!params.overwrite){
    // Check to see if each of the output dirs already exist
    def path_to_check = new File(params.output_dir);
    if (path_to_check.exists()) {
        throw new Exception("""The output directory ${path_to_check} already exists.\n
        If you want to overwrite the contents of this directory, please provide the --overwrite flag.\n
        Alternatively provide --iteration <int> and the output directories will automatically\n
        be suffixed with '_<int>'.""")
    }
}

// Create a channel that is a tuple of the pair_id and the corresponding bam as output from the preprocessing workflow
// The deduction of the pair_ids here is reliant on the bam files ending in ".merged.deduplicated.sorted.bam".
// This will be the case if these bam files have been output from the pre_processing pipeline.
// Alternatively the command line parameter bam_common_extension can be set to a different default string.
Channel.fromFilePairs("${params.bam_input_dir}/*.bam{,.bai}").toList().flatMap{
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
    }.into{gatk_haplotype_caller_gvcf_ch; make_bqsr_tables_bam_ch}

// This scaffold list is created from the reference fasta and is used for the scatter-gather approach
// used in vcf calling.
scaffold_list = {  
    def scaffolds = []
    new File(params.ref_assembly_path).eachLine {
        line -> 
        if (line.startsWith(">")){scaffolds << line.split()[0][1..-1];}
        }
    return scaffolds
}()

// START OF GATK HAPLOTYPE CALLING
process index_dictionary_refgenome{
    container 'broadinstitute/gatk:4.2.0.0'

    input:
    path ref_genome from params.ref_assembly_path

    output:
    tuple path(ref_genome), path("*.dict"), path("*.fai") into gatk_haplotype_caller_ref_genome_ch,GenotypeGVCFs_ref_genome_ch,bcftools_vcfstats_ref_genome_ch

    script:
    """
    gatk CreateSequenceDictionary -R ${ref_genome}
    samtools faidx ${ref_genome}
    """
}

// // Call variants on a per scaffold (chromosone) basis and then gather back to a single vcf (scatter-gather)
// // 1 - Call haplotypes on a per sample basis
// // 2 - run GenomicsDBImport on a per scaffold basis
// // 3 - run GenotypeGVCFs on a per scaffold basis.
// // 4 - Gather the per scaffold files using GatherVcfs
// // 5 - Evaluate the vcfs using bcftools stats and rtg vcfstats
// NB Good reference for the scatter-gather strategy employed here: https://github.com/IARCbioinfo/gatk4-GenotypeGVCFs-nf/blob/master/gatk4-GenotypeGVCFs.nf
// NB HaplotypCaller requires a .fai
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
process gatk_haplotype_caller_gvcf{
    tag {pair_id}
    container 'broadinstitute/gatk:4.2.0.0'
    cpus params.gatk_haplotype_caller_cpus
    memory "24 GB"

    input:
    tuple val(pair_id), path(merged), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from gatk_haplotype_caller_gvcf_ch.combine(gatk_haplotype_caller_ref_genome_ch)

    output:
    tuple val(pair_id), path("${pair_id}.merged.g.vcf.gz") into genomics_db_import_ch
    tuple val(pair_id), path("${pair_id}.merged.g.vcf.gz.tbi") into genomics_db_import_tbi_ch

    script:
    """
    gatk HaplotypeCaller --native-pair-hmm-threads ${task.cpus} -R $ref_genome -I ${merged[0]} -O ${pair_id}.merged.g.vcf.gz -ERC GVCF
    """
}

process genomics_db_import{
    tag "${scaffold}"
    cpus 1
    memory "24 GB"
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

// NB an alternative to passing int the ref_genome_fai and ref_genome_dict
// is to pass the param to the actual path
process GenotypeGVCFs{
    container 'broadinstitute/gatk:4.2.0.0'
	cpus 5
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
    publishDir gatk_output_vcf_publishDir, mode: 'copy'

    input:
    val(scaffolds) from Channel.fromList(scaffold_list).collect()
    path(vcf) from gather_vcfs_for_eval_ch.collect()
	path(vcf_idx) from gather_vcfs_idx_for_eval_ch.collect()

	output:
    tuple path("fowl.eval.vcf"), path("fowl.eval.vcf.idx") into bcftools_stats_ch,rtg_vcfstats_ch

    script:
	"""
	gatk GatherVcfs ${scaffolds.collect{ "--INPUT GenotypeGVCFs.out.${it}.vcf " }.join()} --OUTPUT fowl.eval.vcf
	"""
}

// useful reference for rtg: https://genomics.sschmeier.com/ngs-variantcalling/index.html
process bcftools_vcfstats{
    tag "bcftools_stats"
    container "halllab/bcftools:v1.9"
    publishDir gatk_vcf_stats_publishDir

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
// It outputs a file rtg.stats.summary*.txt (with optional iteration annotation) that contains the total
// SNPs and weighted averges (SNPs per sample) of a selection of metrics including Ti/Tv ratio and insertion/deletion ratio
process rtg_vcfstats_per_sample{
    tag "rtg_vcfstats"
    container 'realtimegenomics/rtg-tools:latest'
    publishDir gatk_vcf_stats_publishDir

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
    publishDir gatk_vcf_stats_publishDir

    input:
    path(rtg_vcfstats_output) from rtg_vcfstats_out_ch

    output:
    path("rtg.stats.summary*.txt") into rtg_vcfstats_summary_out_ch

    script:
    """
    python3 ${bin_dir}/summarise_rtg_vcfstats.py $rtg_vcfstats_output
    """
}