#!/usr/bin/env nextflow

/*
This workflow starts either from a set of .bam/.bam.bai files or per scaffold vcf files.
The input will depend on when this pipeline is being run. If being run directly after producing a set of .bam/.bam/bai files
then the former start point should be used by specifying --input_bam_directory <directory containing .bam/.bam.bai files> on the command line.

If this pipeline is being run after the gatk variant calling pipeline,
then --scaffold_vcf_input_dir <directory containing the .vcf/.vcf.idx files> should be passed at the command line.
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

// Call variants on a per scaffold (chromosone) basis.
// Follow a scatter-gather strategy 
// (first two steps can be skipped if provided with the per scaffold vcfs as produced by the gatk variant calling pipeline)
// 1 - run GenomicsDBImport on a per scaffold basis
// 2 - run GenotypeGVCFs on a per scaffold basis.
// 3 - Hard filter the vcfs on a per scaffold basis using VariantFiltration and remove the filtered variants using SelectVariants
// 4 - Gather the per scaffold files using GatherVcfs GatherVcfs to collect the per chromosome files into a single vcf
// 5 - Make the BQSR recallibration table on a persample basis with BaseRecalibrator using the vcf from 5 as known variants
// 6 - Apply the recalibration table from 6 on a per sample basis with  ApplyBQSR
// 7 - Make a second recalibration table, using the new bam files from 7 as input and the output from 5 as known variants using BaseRecalibrator
// 8 - Produce metrics for the change between the bam files using before and after recalibration using the two tables as input using AnalyzeCovariates
// Good reference for the nextflow coding: https://github.com/IARCbioinfo/gatk4-GenotypeGVCFs-nf/blob/master/gatk4-GenotypeGVCFs.nf

process index_dictionary_refgenome{
    container 'broadinstitute/gatk:latest'

    input:
    path ref_genome from params.ref_assembly_path

    output:
    tuple path(ref_genome), path("*.dict"), path("*.fai") into gatk_haplotype_caller_ref_genome_ch,GenotypeGVCFs_ref_genome_ch,make_bqsr_tables_ref_genome_ch,apply_bqsr_tables_ref_genome_ch,make_bqsr_tables_ref_genome_round_2_ch

    script:
    """
    gatk CreateSequenceDictionary -R ${ref_genome}
    samtools faidx ${ref_genome}
    """
    }

params.scaffold_vcf_input_dir = false
if (!params.scaffold_vcf_input_dir){
    // Then we are working from bams only and need to run all processes.
    // Create a channel that is a tuple of the pair_id and the corresponding bam as output from the preprocessing workflow
    // The deduction of the pair_ids here is reliant on the bam files ending in ".merged.deduplicated.sorted.bam".
    // This will be the case if these bam files have been output from the pre_processing pipeline.
    // Alternatively the command line parameter bam_common_extension can be set to a different default string.
    Channel.fromFilePairs("${params.input_bam_directory}/*.bam{,.bai}").map{it -> [it[1][0].getName().replaceAll(params.bam_common_extension, ""), [it[1][0], it[1][1]]]}.into{gatk_haplotype_caller_gvcf_ch; make_bqsr_tables_bam_ch}
    
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
        gatk HaplotypeCaller --native-pair-hmm-threads ${task.cpus} -R $ref_genome -I ${merged[0]} -O ${pair_id}.merged.g.vcf.gz -ERC GVCF
        """
    }

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

    process GenotypeGVCFs{
        container 'broadinstitute/gatk:latest'
        cpus 5
        tag "${scaffold}"
        publishDir gatk_per_scaffold_vcf_publishDir, mode: 'copy', pattern: '*.{vcf,idx}'

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

}else{
    // Then we have been provided with the per scaffold vcfs as well as the bams so can skip HaplotypeCaller, GenomicsDBImport and GenotypeGVCFs
    // We want to read them into tuples of the scaffold name, the vcf and the vcf.idx
    make_bqsr_tables_bam_ch = Channel.fromFilePairs("${params.input_bam_directory}/*.bam{,.bai}").map{it -> [it[1][0].getName().replaceAll(params.bam_common_extension, ""), [it[1][0], it[1][1]]]}
    hard_filter_ch = Channel.fromFilePairs("${params.scaffold_vcf_input_dir}/*.vcf{,.idx").map{it -> it.flatten()}
}

process hard_filter{
	tag {scaffold}
    container 'broadinstitute/gatk:latest'
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

// NB GatherVcfs must be provided with the vcf files in order of the scaffolds. We do this using the scaffhold_list.
process gather_vcfs{
	tag "GatherVcfs"
    container 'broadinstitute/gatk:latest'

    input:
    val(scaffhold_list) from scaffold_list_gather_vcfs_ch.collect()
    path(vcf) from vcf_gather_vcfs_ch.collect()
	path(vcf_idx) from vcf_idx_gather_vcfs_ch.collect()

	output:
    tuple path("fowl.filtered.vcf"), path("fowl.filtered.vcf.idx") into make_bqsr_tables_known_variants_ch,make_bqsr_tables_known_variants_round_2_ch

    script:
	"""
	gatk GatherVcfs ${scaffhold_list.collect{ "--INPUT ${it}.filtered.vcf " }.join()} --OUTPUT fowl.filtered.vcf
	"""
}


// We will want to do this on a per sample basis
// We will want to use the fowl.vcf and the fowl.vcf.idx multiple times withhout using them up
// It may be that we have to output these files as values rather than path/files in the gather_vcf process.
// make_bqsr_tables_known_variants_ch.view()
process make_bqsr_tables{
    tag {pair_id}
    container 'broadinstitute/gatk:latest'

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
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(pair_id), path(merged_bam), path(table), path(ref_genome), path(ref_genome_dict), path(ref_genome_fai) from apply_bqsr_tables_ch.combine(apply_bqsr_tables_ref_genome_ch)

    output:
    tuple val(pair_id), path("${pair_id}.recalibrated.bam{,.bai}") into apply_bqsr_tables_out_ch,apply_bqsr_tables_round_2_ch

    script:
    """
    gatk ApplyBQSR -R $ref_genome -I ${merged_bam[0]} \
    --bqsr-recal-file $table \
    -O ${pair_id}.recalibrated.bam;
    """
}

// See the below issue to get an idea of the general strategy of getting AnalyzeCovariates to work in 
// GATK 4.
// https://github.com/broadinstitute/gatk/issues/322

process make_bqsr_tables_round_2{
    tag {pair_id}
    container 'broadinstitute/gatk:latest'

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
    container 'broadinstitute/gatk:latest'
    publishDir analyze_covariates_publishDir

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