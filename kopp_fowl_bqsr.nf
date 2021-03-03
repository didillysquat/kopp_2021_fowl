#!/usr/bin/env nextflow

/*
This workflow starts either from a set of .bam/.bam.bai files or per scaffold vcf files.
The input will depend on when this pipeline is being run. If being run directly after producing a set of .bam/.bam/bai files
then the former start point should be used by specifying --bam_input <directory containing .bam/.bam.bai files> on the command line.

If this pipeline is being run after the gatk variant calling pipeline,
then --scaffold_vcf_input <directory containing the .vcf/.vcf.idx files> should be passed at the command line.
This pipeline takes these per scaffold vcf files as input to save having to recompute the
GenomicsDBImport and GenotypeGVCFs processess that also have to be called as part of the variant calling pipeline.

This pipeline outputs a set of .bam/.bam.bai files that have been through BQSR.
It performs 1 round of BQSR everytime it is run.
Every time it runs it produces a set of metrics that charactersise the difference between
the before and after BQSR bam sets.
*/

if (params.subsample){
    params.output_dir = "${workflow.launchDir}/outputs_sub_sampled/"
    gatk_per_scaffold_vcf_publishDir = [params.output_dir, "gvcfs_sub_${params.subsample_depth}"].join(File.separator)
    gatk_output_vcf_publishDir = [params.output_dir, "gatk_output_variants_sub_${params.subsample_depth}"].join(File.separator)
    gatk_vcf_stats_publishDir = [params.output_dir, "vcf_stats_sub_${params.subsample_depth}"].join(File.separator)
}else{
    params.output_dir = "${workflow.launchDir}/outputs"
    gatk_per_scaffold_vcf_publishDir = [params.output_dir, "gvcfs_GVCFs"].join(File.separator)
    gatk_output_vcf_publishDir = [params.output_dir, "gatk_output_variants"].join(File.separator)
    gatk_vcf_stats_publishDir = [params.output_dir, "vcf_stats"].join(File.separator)
}


// Create a channel that is a tuple of the pair_id and the corresponding bam as output from the preprocessing workflow
// The deduction of the pair_ids here is reliant on the bam files ending in ".merged.deduplicated.sorted.bam".
// This will be the case if these bam files have been output from the pre_processing pipeline.
// Alternatively the command line parameter bam_common_extension can be set to a different default string.
Channel.fromFilePairs("${params.input_bam_directory}/*.bam{,*.bai}").map{it -> [it[1][0].getName().replaceAll(params.bam_common_extension, ""), [it[1][0], it[1][1]]]}.into{gatk_haplotype_caller_gvcf_ch; make_bqsr_tables_bam_ch}

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
    container 'broadinstitute/gatk:latest'

    input:
    path ref_genome from params.ref_assembly_path

    output:
    tuple path(ref_genome), path("*.dict"), path("*.fai") into gatk_haplotype_caller_ref_genome_ch

    script:
    """
    gatk CreateSequenceDictionary -R ${ref_genome}
    samtools faidx ${ref_genome}
    """
}

// NB HaplotypCaller requires a .fai
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
process gatk_haplotype_caller_gvcf{
    tag {pair_id}
    container 'broadinstitute/gatk:latest'
    cpus params.gatk_haplotype_caller_cpus
    memory "24 GB"

    input:
    tuple val(pair_id), path(merged), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from gatk_haplotype_caller_gvcf_ch.combine(gatk_haplotype_caller_ref_genome_ch)

    output:
    tuple val(pair_id), path("${pair_id}.merged.g.vcf.gz") into genomics_db_import_ch
    tuple val(pair_id), path("${pair_id}.merged.g.vcf.gz.tbi") into genomics_db_import_tbi_ch

    script:
    """
    gatk HaplotypeCaller --native-pair-hmm-threads ${task.cpus} -R ${params.ref_assembly_path} -I ${merged[0]} -O ${pair_id}.merged.g.vcf.gz -ERC GVCF
    """
}

// Call variants on a per scaffold (chromosone) basis.
// Follow a scatter-gather strategy 
// 1 - run GenomicsDBImport on a per scaffold basis
// 2 - run GenotypeGVCFs on a per scaffold basis.
// 3 - Hard filter the vcfs on a per scaffold basis using VariantFiltration and remove the filtered variants using SelectVariants
// 4 - Gather the per scaffold files using GatherVcfs GatherVcfs to collect the per chromosome files into a single vcf
// 5 - Make the BQSR recallibration table on a persample basis with BaseRecalibrator using the vcf from 5 as known variants
// 6 - Apply the recalibration table from 6 on a per sample basis with  ApplyBQSR
// 7 - Make a second recalibration table, using the new bam files from 7 as input and the output from 5 as known variants using BaseRecalibrator
// 8 - Produce metrics for the change between the bam files using before and after recalibration using the two tables as input using AnalyzeCovariates
// Good reference for the nextflow coding: https://github.com/IARCbioinfo/gatk4-GenotypeGVCFs-nf/blob/master/gatk4-GenotypeGVCFs.nf
process genomics_db_import{
    tag "${scaffold}"
    cpus 1
    memory "24 GB"
    container 'broadinstitute/gatk:latest'

    input:
    each scaffold from Channel.fromList(scaffold_list)
    path(gvcf) from genomics_db_import_ch.collect{it[1]}
    path(gvcf_tbi) from genomics_db_import_tbi_ch.collect{it[1]}
	
	output:
    tuple val(scaffold), file("genomicsdbi.out.${scaffold}") into genotype_GVCFs_ch
	
    script:
	"""
	gatk  --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenomicsDBImport ${gvcf.collect { "-V $it " }.join()} --genomicsdb-workspace-path genomicsdbi.out.${scaffold} -L $scaffold --batch-size 50
	"""
}

// NB an alternative to passing int the ref_genome_fai and ref_genome_dict
// is to pass the param to the actual path
process GenotypeGVCFs{
    container 'broadinstitute/gatk:latest'
	cpus 5
	tag "${scaffold}"
	publishDir gatk_per_scaffold_vcf_publishDir, mode: 'copy', pattern: '*.{vcf,idx}'

    input:
	tuple val(scaffold), file(workspace) from genotype_GVCFs_ch
   	path genome from params.ref_assembly_path
    path ref_genome_fai from ref_assembly_fai_path
    path ref_genome_dict from ref_assembly_dict_path

	output:
    tuple val(scaffold), file("GenotypeGVCFs.out.${scaffold}.vcf"), file("GenotypeGVCFs.out.${scaffold}.vcf.idx") into hard_filter_ch
    file("GenotypeGVCFs.out.${scaffold}.vcf") into gather_vcfs_for_eval_ch
    file("GenotypeGVCFs.out.${scaffold}.vcf.idx") into gather_vcfs_idx_for_eval_ch

    script:
	"""
    WORKSPACE=\$( basename ${workspace} )
    gatk GenotypeGVCFs -R ${genome} -O GenotypeGVCFs.out.${scaffold}.vcf \
    --only-output-calls-starting-in-intervals -V gendb://\$WORKSPACE -L ${scaffold}
	"""
}