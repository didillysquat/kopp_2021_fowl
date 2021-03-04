#!/usr/bin/env nextflow

/*
This workflow takes bam files as input and outputs a) multi-sample per scaffold vcfs and b) multi-sample vcfs.
The multi-sample per-scaffold vcfs (a), can be used as input to BQSR workflow to save running GenomicsDBImport and GenotypeGVCFs
as the first two processes of the BQSR pipeline.
The multi-sample vcfs (b) are the vcfs that will be carried forwards for analysis.
This workflow outputs summary metrics of the produced multi-sample vcf files (b) that can be used to compare the effect
running BQSR (e.g. by comparing Ti/Tv ratios or number of SNPs).
*/

bin_dir = "${workflow.launchDir}/bin"

if (params.subsample){
    params.output_dir = "${workflow.launchDir}/outputs_sub_sampled/"
    gatk_per_scaffold_vcf_publishDir = [params.output_dir, "gatk_scaffold_vcfs_sub_${params.subsample_depth}"].join(File.separator)
    gatk_output_vcf_publishDir = [params.output_dir, "gatk_output_variants_sub_${params.subsample_depth}"].join(File.separator)
    gatk_vcf_stats_publishDir = [params.output_dir, "gatk_vcf_stats_sub_${params.subsample_depth}"].join(File.separator)
}else{
    params.output_dir = "${workflow.launchDir}/outputs"
    gatk_per_scaffold_vcf_publishDir = [params.output_dir, "gatk_scaffold_vcfs_GVCFs"].join(File.separator)
    gatk_output_vcf_publishDir = [params.output_dir, "gatk_output_variants"].join(File.separator)
    gatk_vcf_stats_publishDir = [params.output_dir, "gatk_vcf_stats"].join(File.separator)
}

// Create a channel that is a tuple of the pair_id and the corresponding bam as output from the preprocessing workflow
// The deduction of the pair_ids here is reliant on the bam files ending in ".merged.deduplicated.sorted.bam".
// This will be the case if these bam files have been output from the pre_processing pipeline.
// Alternatively the command line parameter bam_common_extension can be set to a different default string.
Channel.fromFilePairs("${params.input_bam_directory}/*.bam{,.bai}").map{it -> [it[1][0].getName().replaceAll(params.bam_common_extension, ""), [it[1][0], it[1][1]]]}.into{gatk_haplotype_caller_gvcf_ch; make_bqsr_tables_bam_ch}

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
    tuple path(ref_genome), path("*.dict"), path("*.fai") into gatk_haplotype_caller_ref_genome_ch,GenotypeGVCFs_ref_genome_ch,bcftools_vcfstats_ref_genome_ch

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
    gatk HaplotypeCaller --native-pair-hmm-threads ${task.cpus} -R $ref_genome -I ${merged[0]} -O ${pair_id}.merged.g.vcf.gz -ERC GVCF
    """
}

// // Call variants on a per scaffold (chromosone) basis and then gather back to a single vcf (scatter-gather)
// // 1 - run GenomicsDBImport on a per scaffold basis
// // 2 - run GenotypeGVCFs on a per scaffold basis.
// // 3 - Gather the per scaffold files using GatherVcfs
// // 4 - Evaluate the vcfs using bcftools stats and rtg vcfstats
// Good reference for the scatter-gather strategy employed here: https://github.com/IARCbioinfo/gatk4-GenotypeGVCFs-nf/blob/master/gatk4-GenotypeGVCFs.nf
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
    tuple val(scaffold), path("genomicsdbi.out.${scaffold}") into genotype_GVCFs_ch
	
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
    container 'broadinstitute/gatk:latest'
    publishDir gatk_output_vcf_publishDir, mode: 'copy'

    input:
    val(scaffhold_list) from Channel.fromList(scaffold_list).collect()
    path(vcf) from gather_vcfs_for_eval_ch.collect()
	path(vcf_idx) from gather_vcfs_idx_for_eval_ch.collect()

	output:
    tuple path("fowl.eval.vcf"), path("fowl.eval.vcf.idx") into bcftools_stats_ch,rtg_vcfstats_ch

    script:
	"""
	gatk GatherVcfs ${scaffhold_list.collect{ "--INPUT GenotypeGVCFs.out.${it}.vcf " }.join()} --OUTPUT fowl.eval.vcf
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
    container 'broadinstitute/gatk:latest'
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