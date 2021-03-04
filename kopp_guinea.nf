#!/usr/bin/env nextflow

params.raw_reads_dir = "/home/humebc/projects/20210125_kopp_guinea_fowl/raw_seq_files"
params.ref_assembly_path = "/home/humebc/projects/20210125_kopp_guinea_fowl/nf_pipeline_kopp/nummel_ref_assembly/GCF_002078875.1_NumMel1.0_genomic.fna"
// NB both the .fai file and the .dict file are created in the process index_dictionary_refgenome
// We pass them into the process GenotypeGVCFs as values created directly from these variables so that they can
// be reused rather than pass in the output files from the index_dictionary_refgenome which may cause them to be used up in
// a single process (we have not tested whether this is true but it most likely is).
ref_assembly_fai_path = params.ref_assembly_path+".fai"
ref_assembly_dict_path = params.ref_assembly_path.replaceAll(".fna", ".dict")
File ref_assembly_file = new File(params.ref_assembly_path);
ref_assembly_dir = ref_assembly_file.getParent();

// Check that the ref genome is decompressed
if (params.ref_assembly_path.endsWith(".bgz") || params.ref_assembly_path.endsWith(".gz") || params.ref_assembly_path.endsWith(".zip")){
    throw new Exception("The reference assembly genome must be decompressed and in fasta format (.fna, .fa, .fasta). Currently the reference genome path is set to ${params.ref_assembly_path}")
}

params.path_to_read_group_lists = "/home/humebc/projects/20210125_kopp_guinea_fowl/raw_seq_files/file-list-run-1,/home/humebc/projects/20210125_kopp_guinea_fowl/raw_seq_files/file-list-run-2"

bin_dir = "${workflow.launchDir}/bin"
envs_dir = "${workflow.launchDir}/envs"
launch_dir = "${workflow.launchDir}"
tru_seq_pe_fasta_path = "${workflow.launchDir}/TruSeq3-PE.fa"

params.subsample = true
params.subsample_depth = 10000
// We originally made a map to use in the processes but this ran into 
// multiprocessing conflict issues and made the dictionary unusable unless running
// with maxforks set to 1.
// To be more nextflow, we will make a list of tuplets with pair ID to readgroup string.
// Then we will do a join operation using the pair_id as the keymake 
// exmple string: RGID=xxx RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=xxx
// NB it is critical to convert the key for a map to a string using .toString() if using e.g. "${my_var}"

scaffold_list = {  
    def scaffolds = []
    new File(params.ref_assembly_path).eachLine {
        line -> 
        if (line.startsWith(">")){scaffolds << line.split()[0][1..-1];}
        }
    return scaffolds
}()

Channel.fromList(scaffold_list).into{scaffold_list_ch; scaffold_list_gather_vcfs_ch}
// NB the sample number must be unique for all samples irrespective of read group.
read_group_map = { 
    String[] file_lists = params.path_to_read_group_lists.split(",");
    def tup_list = []
    def j = 0
    file_lists.eachWithIndex{ file_list, i -> 
            // i will be the RGID value
            // j will be the RGSM value
        
        def pair_id_list = []
        new File(file_list).eachLine {  
            line ->
            if (line.endsWith('.gz')){
                
                // The this is a seq file and we want to get its pair_id and add this to the map
                def pattern = ~/^(?<pair>.*)_R[1-2].*$/
                def file_name = line.split('/')[-1]
                def matcher = file_name =~ pattern
                if (params.subsample){
                    
                    if (matcher){
                        pair_id = matcher.group("pair")
                        if (params.subsample){
                            pair_id = "${pair_id}_sub_${params.subsample_depth}"
                        }
                        // At this point we need to see if the mapping information has already been added
                        // If it has not, then add in the mapping information 
                        if (!pair_id_list.contains(pair_id)){
                            pair_id_list << pair_id
                            j++
                            tup_list << ["${pair_id}".toString(), "RGID=${i+1} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${j}"]
                        }                            
                    }else{
                        throw new Exception("An error has occured when making the readgroup dictionary\nNo match was found for ${file_name}")
                    }
                }
            }
        }
    }
    return tup_list     
}()

// Publish dirs
if (params.subsample){
    params.output_dir = "${workflow.launchDir}/outputs_sub_sampled/"
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim_sub_${params.subsample_depth}"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim_sub_${params.subsample_depth}"].join(File.separator)
    markduplicates_metrics_publishDir = [params.output_dir, "markduplicates_metrics_sub_${params.subsample_depth}"].join(File.separator)
    pre_seq_c_curve_publish_dir = [params.output_dir, "pre_seq_c_curve_sub_${params.subsample_depth}"].join(File.separator)
    collect_gc_bias_metrics_publishDir = [params.output_dir, "collect_gc_bias_metrics_sub_${params.subsample_depth}"].join(File.separator)
    pcr_bottleneck_coefficient_publishDir = [params.output_dir, "pcr_bottleneck_coefficient_sub_${params.subsample_depth}"].join(File.separator)
    mosdepth_sequencing_coverage_publishDir = [params.output_dir, "mosdepth_sequencing_coverage_sub_${params.subsample_depth}"].join(File.separator)
    genotype_GVCFs_publishDir = [params.output_dir, "genotype_GVCFs_sub_${params.subsample_depth}"].join(File.separator)
    gatk_output_variants_publishDir = [params.output_dir, "gatk_output_variants_sub_${params.subsample_depth}"].join(File.separator)
    vcf_stats_publishDir = [params.output_dir, "vcf_stats_sub_${params.subsample_depth}"].join(File.separator)
    analyze_covariates_publishDir = [params.output_dir, "analyze_covariates_sub_${params.subsample_depth}"].join(File.separator)
    
}else{
    params.output_dir = "${workflow.launchDir}/outputs"
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim"].join(File.separator)
    markduplicates_metrics_publishDir = [params.output_dir, "markduplicates_metrics"].join(File.separator)
    pre_seq_c_curve_publish_dir = [params.output_dir, "pre_seq_c_curve"].join(File.separator)
    collect_gc_bias_metrics_publishDir = [params.output_dir, "collect_gc_bias_metrics"].join(File.separator)
    pcr_bottleneck_coefficient_publishDir = [params.output_dir, "pcr_bottleneck_coefficient"].join(File.separator)
    mosdepth_sequencing_coverage_publishDir = [params.output_dir, "mosdepth_sequencing_coverage"].join(File.separator)
    genotype_GVCFs_publishDir = [params.output_dir, "genotype_GVCFs"].join(File.separator)
    gatk_output_variants_publishDir = [params.output_dir, "gatk_output_variants"].join(File.separator)
    vcf_stats_publishDir = [params.output_dir, "vcf_stats"].join(File.separator)
    analyze_covariates_publishDir = [params.output_dir, "analyze_covariates"].join(File.separator)
}

// CPUs
params.trimmomatic_threads = 50
params.nextgenmap_threads = 40
params.markduplicates_cpus = 20
params.gatk_haplotype_caller_cpus = 5

// START OF PRE PROCESSING
/* 
If subsample, create a channel that will pass into the subsampling process and then pass the subsampled
files into the trimming and fastqc
Modify the pair name to indicate subsampled files. This way it will carry through the remainder of the anlysis.
if not subsample, then create channel directly from the raw sequencing files
*/
if (params.subsample){
    ch_subsample = Channel.fromFilePairs("${params.raw_reads_dir}/*R{1,2}*.fastq.gz")
    // ch_ref_genome = Channel.fromPath(params.ref_assembly_path)
    ch_ref_genome = Channel.value(params.ref_assembly_path)
    
    process subsample{
        tag "${pair_id}"
        cache 'lenient'
        container 'biocontainers/seqtk:v1.3-1-deb_cv1'
        containerOptions '-u $(id -u):$(id -g)'
        publishDir "${params.output_dir}/${params.subsample_depth}_subsampled_reads", mode: 'copy'

        input:
        tuple val(pair_id), file(reads) from ch_subsample

        output:
        tuple val("${pair_id}_sub_${params.subsample_depth}"), file("${pair_id}_sub_${params.subsample_depth}*.fastq.gz") into ch_fastqc_pre_trim,ch_trimmomatic_input

        script:
        read_out_one = reads[0].getName().replaceAll("${pair_id}", "${pair_id}_sub_${params.subsample_depth}")
        read_out_two = reads[1].getName().replaceAll("${pair_id}", "${pair_id}_sub_${params.subsample_depth}")
        
        """
        seqtk sample -s100 ${reads[0]} ${params.subsample_depth} | gzip > ${read_out_one}
        seqtk sample -s100 ${reads[1]} ${params.subsample_depth} | gzip > ${read_out_two}
        """
    }
}else{
    Channel.fromFilePairs("${params.raw_reads_dir}/*_{1,2}.fastq.gz").into{ch_fastqc_pre_trim; ch_trimmomatic_input}
    // ch_ref_genome = Channel.fromPath(params.ref_assembly_path)
    // ch_ref_genome_name = Channel.value(params.ref_assembly_path)
}

process fastqc_pre_trim{
    tag "${pair_id}"
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir fastqc_pre_trim_publish_dir, mode: 'copy'

    input:
    tuple val(pair_id), file(fastq_file) from ch_fastqc_pre_trim.flatMap{[["${it[0]}_1", it[1][0]], ["${it[0]}_2", it[1][1]]]}

    output:
    file "${out_name}" into ch_fastqc_pre_trim_output

    script:
    out_name = fastq_file.getName().replaceAll('.fastq.gz', '.pre_trim.fastqc.html')
    """
    fastqc -o . $fastq_file
    mv *.html ${out_name}
    """
}

process trimmomatic{
    cache 'lenient'
	tag "${pair_id}"
    container 'davelabhub/trimmomatic:0.39--1'
    containerOptions '-u $(id -u):$(id -g)'
	
	input:
	tuple val(pair_id), file(fastqs) from ch_trimmomatic_input
	
	output:
	tuple val(pair_id), file("${pair_id}*{1,2}P.fq.gz") into ch_fastqc_post_trim_paired,ch_ngm_paired
    tuple val(pair_id), file("${pair_id}*{1,2}U.fq.gz") into ch_fastqc_post_trim_unpaired,ch_ngm_unpaired

	script:
	outbase = fastqs[0].getName().replaceAll('.fastq.gz', '.trimmed.fq.gz')
	"""
	trimmomatic PE -threads ${params.trimmomatic_threads} -basein ${fastqs[0]} \\
		-baseout $outbase \\
		ILLUMINACLIP:${tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
		SLIDINGWINDOW:4:20 MINLEN:70 LEADING:10 HEADCROP:5 TRAILING:10 
	"""
}


process fastqc_post_trim{
    tag "${pair_id}"
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir fastqc_post_trim_publish_dir, mode: 'copy'

    input:
    tuple val(pair_id), file(paired_files), file(unpaired_files) from ch_fastqc_post_trim_paired.join(ch_fastqc_post_trim_unpaired)

    output:
    tuple file("${pair_id}*1P*fastqc.html"), file("${pair_id}*2P*fastqc.html"), file("${pair_id}*1U*fastqc.html"), file("${pair_id}*2U*fastqc.html") into ch_fastqc_post_trim_output

    script:
    """
    fastqc -o . ${paired_files[0]}
    fastqc -o . ${paired_files[1]}
    fastqc -o . ${unpaired_files[0]}
    fastqc -o . ${unpaired_files[1]}
    """
}

// Index the nextgenmap index before doing the mapping
// If the reference is already indexed, the index will not happen again
// NB We invested a lot of time into getting the flow of this indexing process to work
// (i.e. as a precursor to the mapping). The key is to use a value channel (rather than a queue channel)
// or to directly create the path object from the param variable. This then means that the
// value of this channel can be used multiple times (once for each mapping procedure).
process nextgenmap_indexing{
    conda 'envs/ngm.yaml'
    cpus params.nextgenmap_threads

    input:
    path ref_genome from params.ref_assembly_path

    output:
    path ref_genome into index_ch

    script:
    """
     ngm -t ${task.cpus} -r ${params.ref_assembly_path}
    """
}

// // // We don't use the genome variable that comes from the index_ch
// // // This input is only there to ensure that the indexing has been performed before the mapping
process nextgenmap_mapping{
    tag "${pair_id}"
    conda 'envs/ngm.yaml'
    cpus params.nextgenmap_threads

    input:
    path genome from index_ch
	tuple val(pair_id), file(paired), file(unpaired) from ch_ngm_paired.join(ch_ngm_unpaired)

    output:
    tuple val(pair_id), file("${pair_id}_P_mapped.sam"), file("${pair_id}_U1_mapped.sam"), file("${pair_id}_U2_mapped.sam") into add_read_group_headers_ch

    script:
    """
    ngm -t ${task.cpus} -r ${params.ref_assembly_path} -1 ${paired[0]} -2 ${paired[1]} -o ${pair_id}_P_mapped.sam;
    ngm -t ${task.cpus} -r ${params.ref_assembly_path} -q ${unpaired[0]} -o ${pair_id}_U1_mapped.sam;
    ngm -t ${task.cpus} -r ${params.ref_assembly_path} -q ${unpaired[1]} -o ${pair_id}_U2_mapped.sam;
    """
}


process merge_paired_and_unpaired{
    tag "${pair_id}"
    container 'broadinstitute/gatk:4.1.9.0'

    input:
    tuple val(pair_id), file(paired), file(unpaired_one), file(unpaired_two) from add_read_group_headers_ch

    output:
    tuple val(pair_id), file("${pair_id}_all_reads_mapped.sam") into read_groups_all_ch

    script:
    """
    gatk MergeSamFiles --INPUT $paired --INPUT $unpaired_one --INPUT $unpaired_one --OUTPUT ${pair_id}_all_reads_mapped.sam
    """

}

process add_read_group_headers{
    tag "${pair_id}"
    container 'broadinstitute/gatk:4.1.9.0'

    input:
    tuple val(pair_id), file(merged), val(read_group_string) from read_groups_all_ch.join(Channel.fromList(read_group_map))

    output:
    tuple val(pair_id), file("${pair_id}.merged.mapped.headers.sam") into mark_duplicates_ch

    script:
    """
    gatk AddOrReplaceReadGroups I=${merged} O=${pair_id}.merged.mapped.headers.sam ${read_group_string}
    """
}

process markduplicates_spark{
    tag "${pair_id}"
    container 'broadinstitute/gatk:4.1.9.0'
    publishDir markduplicates_metrics_publishDir, pattern: "*.metrics.txt"

    input:
    tuple val(pair_id), file(merged) from mark_duplicates_ch

    output:
    tuple val(pair_id), file("${pair_id}.merged.deduplicated.sorted.bam{,.bai}") into pre_seq_c_curve_ch,collect_gc_bias_metrics_ch,pcr_bottleneck_coefficient_ch,mosdepth_sequencing_coverage_ch,gatk_haplotype_caller_gvcf_ch,make_bqsr_tables_bam_ch
    file("${pair_id}.merged.deduplicated.sorted.metrics.txt") into mark_duplicate_metrics_ch

    script:
    """
    gatk MarkDuplicatesSpark --create-output-bam-index --remove-sequencing-duplicates -I ${merged} -O ${pair_id}.merged.deduplicated.sorted.bam -M ${pair_id}.merged.deduplicated.sorted.metrics.txt
    """
}
// END OF PRE PROCESSING

// START OF EVALUATION METRICS
process collect_gc_bias_metrics{
    tag pair_id
    publishDir collect_gc_bias_metrics_publishDir, mode: 'copy'
    container 'broadinstitute/gatk:4.1.9.0'

    input:
    tuple val(pair_id), path(merged) from collect_gc_bias_metrics_ch
    path ref_genome from params.ref_assembly_path
    
    output:
    tuple path("${pair_id}.merged.GCBias.txt"), path("${pair_id}.merged.GCBias.pdf"), path("${pair_id}.merged.SumBias.txt")

    script:
    """
    gatk CollectGcBiasMetrics I=${merged[0]} O=${pair_id}.merged.GCBias.txt \
	CHART=${pair_id}.merged.GCBias.pdf S=${pair_id}.merged.SumBias.txt R=${ref_genome}
    """
}


process pcr_bottleneck_coefficient{
    tag pair_id
    publishDir pcr_bottleneck_coefficient_publishDir, mode: 'copy'
    container 'encodedcc/atac-seq-pipeline:PIP-1469_pbam_b92239a0-82a0-4297-b8b3-3a655a4626a8'

    input:
    tuple val(pair_id), path(merged) from pcr_bottleneck_coefficient_ch
    
    output:
    tuple path("${pair_id}.merged.allout.txt"), path("${pair_id}.pbc.merged.txt") into pcr_bottleneck_coefficient_out_ch

    shell:
    '''
    bedtools bamtobed -bedpe -i !{merged[0]} | awk 'BEGIN{{OFS="\\t"}}{{print $1, $2, $4, $6, $9, $10}}' | grep -v "^{}\\s" | sort | uniq -c  | \
    awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1) {{m1=m1+1}} ($1==2) {{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{m1_m2=-1.0;
    if (m2>0) m1_m2=m1/m2;
    m0_mt=0;
    if (mt>0) m0_mt=m0/mt;
    m1_m0=0;
    if (m0>0) m1_m0=m1/m0;
    printf "%d %d %d %d %f %f %f\\n",mt,m0,m1,m2,m0_mt,m1_m0,m1_m2}}\' > l.txt;
    cat l.txt >> !{pair_id}.merged.allout.txt;
    awk '{print $5}' l.txt > a.txt;
    cat a.txt >> !{pair_id}.pbc.merged.txt;
    '''
}

// -n = dont output per-base depth. skipping this output will speed execution
// -x = dont look at internal cigar operations or correct mate overlaps (recommended for most use-cases)
process mosdepth_sequencing_coverage{
    tag pair_id
    publishDir mosdepth_sequencing_coverage_publishDir, mode: 'copy'
    container 'davelabhub/mosdepth:0.2.5--hb763d49_0'
    
    input:
    tuple val(pair_id), path(merged) from mosdepth_sequencing_coverage_ch
    
    output:
    tuple val(pair_id), path("${pair_id}.merged.mosdepth.global.dist.txt") into mosdepth_sequencing_coverage_out_ch

    script:
    """
    mosdepth -nx ${pair_id}.merged ${merged[0]}
    """
}


process mosdepth_plot_seq_coverage{
    tag pair_id
    publishDir mosdepth_sequencing_coverage_publishDir, mode: 'copy'
    conda "${envs_dir}/python3.yaml"

    input:
    tuple val(pair_id), path(merged) from mosdepth_sequencing_coverage_out_ch

    output:
    path("${pair_id}.merged.mostdepth.html") into mosdepth_plot_seq_coverage_out_ch

    script:
    """
    python3 ${bin_dir}/plot-dist.py -o ${pair_id}.merged.mostdepth.html ${merged}
    """
}
// END OF EVALUATION METRICS

// START OF GATK HAPLOTYPE CALLING
// // TODO the copy can be changed to move at a later date.
// // If we use the ref_genome as the path to the ref genome, then the output is in the working dir of this
// // nextflow process. But if we use params.ref_assembly_path then the .fai and .dict are output in
// // the params.ref_asssembly_path directory. WHich is good for now but could cause complicaitons
// // once were on BINAC or some
// process index_dictionary_refgenome{
//     container 'broadinstitute/gatk:4.1.9.0'
//     publishDir ref_assembly_dir, mode: 'copy', saveAs: {filename -> if(filename.toString().endsWith(".fai") || filename.toString().endsWith(".dict")){return "${filename}";}else{return null;}}

//     input:
//     path ref_genome from params.ref_assembly_path

//     output:
//     tuple path("*.dict"), path("*.fai") into index_dictionary_refgenome_out_ch
//     path ref_genome into gatk_index_ch

//     script:
//     """
//     gatk CreateSequenceDictionary -R ${ref_genome}
//     samtools faidx ${ref_genome}
//     """
// }


// // // NB HaplotypCaller requires a .fai
// // // https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
// process gatk_haplotype_caller_gvcf{
//     tag pair_id
//     container 'broadinstitute/gatk:4.1.9.0'
//     cpus params.gatk_haplotype_caller_cpus

//     input:
//     tuple val(pair_id), path(merged) from gatk_haplotype_caller_gvcf_ch
//     path ref_genome from gatk_index_ch

//     output:
//     tuple val(pair_id), path("${pair_id}.merged.g.vcf.gz") into genomics_db_import_ch
//     tuple val(pair_id), path("${pair_id}.merged.g.vcf.gz.tbi") into genomics_db_import_tbi_ch

//     script:
//     """
//     gatk HaplotypeCaller --native-pair-hmm-threads ${task.cpus} -R ${params.ref_assembly_path} -I ${merged[0]} -O ${pair_id}.merged.g.vcf.gz -ERC GVCF
//     """
// }

// // // Call variants on a per scaffold (chromosone) basis.
// // // Follow a scatter-gather strategy 
// // // 1 - run GenomicsDBImport on a per scaffold basis
// // // 2 - run GenotypeGVCFs on a per scaffold basis.
// // // 3 - Gather the per scaffold files using GatherVcfs and evaluate using bcftools stats and rtg vcfstats
// // // 4 - Hard filter the vcfs on a per scaffold basis using VariantFiltration and remove the filtered variants using SelectVariants
// // // 5 - Gather the per scaffold files using GatherVcfs GatherVcfs to collect the per chromosome files into a single vcf
// // // 6 - Make the BQSR recallibration table on a persample basis with BaseRecalibrator using the vcf from 5 as known variants
// // // 7 - Apply the recalibration table from 6 on a per sample basis with  ApplyBQSR
// // // 8 - Make a second recalibration table, using the new bam files from 7 as input and the output from 5 as known variants using BaseRecalibrator
// // // 9 - Produce metrics for the change between the bam files using before and after recalibration using the two tables as input using AnalyzeCovariates
// // Good reference for the nextflow coding: https://github.com/IARCbioinfo/gatk4-GenotypeGVCFs-nf/blob/master/gatk4-GenotypeGVCFs.nf
// process genomics_db_import{
//     tag "${scaffold}"
//     cpus 5
//     container 'broadinstitute/gatk:4.1.9.0'

//     input:
//     each scaffold from scaffold_list_ch
//     path(gvcf) from genomics_db_import_ch.collect{it[1]}
//     path(gvcf_tbi) from genomics_db_import_tbi_ch.collect{it[1]}
	
// 	output:
//     tuple val(scaffold), file("genomicsdbi.out.${scaffold}") into genotype_GVCFs_ch
	
//     script:
// 	"""
// 	gatk  --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' GenomicsDBImport ${gvcf.collect { "-V $it " }.join()} --genomicsdb-workspace-path genomicsdbi.out.${scaffold} -L $scaffold --batch-size 50
// 	"""
// }

// // NB an alternative to passing int the ref_genome_fai and ref_genome_dict
// // is to pass the param to the actual path
// process GenotypeGVCFs{
//     container 'broadinstitute/gatk:4.1.9.0'
// 	cpus 5
// 	tag "${scaffold}"

// 	publishDir genotype_GVCFs_publishDir, mode: 'copy', pattern: '*.{vcf,idx}'

//     input:
// 	tuple val(scaffold), file(workspace) from genotype_GVCFs_ch
//    	path genome from params.ref_assembly_path
//     path ref_genome_fai from ref_assembly_fai_path
//     path ref_genome_dict from ref_assembly_dict_path

// 	output:
//     tuple val(scaffold), file("GenotypeGVCFs.out.${scaffold}.vcf"), file("GenotypeGVCFs.out.${scaffold}.vcf.idx") into hard_filter_ch
//     file("GenotypeGVCFs.out.${scaffold}.vcf") into gather_vcfs_for_eval_ch
//     file("GenotypeGVCFs.out.${scaffold}.vcf.idx") into gather_vcfs_idx_for_eval_ch

//     script:
// 	"""
//     WORKSPACE=\$( basename ${workspace} )
//     gatk GenotypeGVCFs -R ${genome} -O GenotypeGVCFs.out.${scaffold}.vcf \
//     --only-output-calls-starting-in-intervals -V gendb://\$WORKSPACE -L ${scaffold}
// 	"""
// }

// process gather_vcfs_for_eval{
// 	tag "GatherVcfs_for_eval"
//     container 'broadinstitute/gatk:4.1.9.0'
//     publishDir gatk_output_variants_publishDir

//     input:
//     val(scaffhold_list) from Channel.fromList(scaffold_list).collect()
//     path(vcf) from gather_vcfs_for_eval_ch.collect()
// 	path(vcf_idx) from gather_vcfs_idx_for_eval_ch.collect()

// 	output:
//     tuple path("fowl.eval.vcf"), path("fowl.eval.vcf.idx") into bcftools_stats_ch,rtg_vcfstats_ch

//     script:
// 	"""
// 	gatk GatherVcfs ${scaffhold_list.collect{ "--INPUT GenotypeGVCFs.out.${it}.vcf " }.join()} --OUTPUT fowl.eval.vcf
// 	"""
// }

// // useful reference for rtg: https://genomics.sschmeier.com/ngs-variantcalling/index.html
// process bcftools_vcfstats{
//     tag "bcftools_stats"
//     container "halllab/bcftools:v1.9"
//     publishDir vcf_stats_publishDir

//     input:
//     tuple path(vcf), path(vcf_idx) from bcftools_stats_ch

//     output:
//     path("bcftools.stats.txt") into bcftools_stats_out_ch

//     script:
//     """
//     bcftools stats -F ${params.ref_assembly_path} -s - $vcf > bcftools.stats.txt
//     """
// }


// // This produces a really handy per sample output but ideally we'd also like a summary metric
// // We have written and use summarise_rtg_vcfstats.py to do this.
// // It outputs a file rtg.stats.summary*.txt (with optional iteration annotation) that contains the total
// // SNPs and weighted averges (SNPs per sample) of a selection of metrics including Ti/Tv ratio and insertion/deletion ratio
// process rtg_vcfstats_per_sample{
//     tag "rtg_vcfstats"
//     container 'realtimegenomics/rtg-tools:latest'
//     publishDir vcf_stats_publishDir

//     input:
//     tuple path(vcf), path(vcf_idx) from rtg_vcfstats_ch

//     output:
//     path("rtg.stats.per_sample.txt") into rtg_vcfstats_out_ch

//     script:
//     """
//     rtg vcfstats $vcf > rtg.stats.per_sample.txt
//     """

// }

// process rtg_vcfstats_summary{
//     tag "rtg_vcfstats_summary"
//     container 'broadinstitute/gatk:4.1.9.0'
//     publishDir vcf_stats_publishDir

//     input:
//     path(rtg_vcfstats_output) from rtg_vcfstats_out_ch

//     output:
//     path("rtg.stats.summary*.txt") into rtg_vcfstats_summary_out_ch

//     script:
//     """
//     python3 ${bin_dir}/summarise_rtg_vcfstats.py $rtg_vcfstats_output
//     """
// }

// process hard_filter{
// 	tag {scaffold}
//     container 'broadinstitute/gatk:4.1.9.0'
//     cpus 5

//     input:
// 	tuple val(scaffold), path(vcf), path(vcfidx) from hard_filter_ch

// 	output:
//     file("${scaffold}.filtered.vcf") into vcf_gather_vcfs_ch
//     file("${scaffold}.filtered.vcf.idx") into vcf_idx_gather_vcfs_ch

//     script:
// 	"""
// 	gatk VariantFiltration \
//     -filter "Qual >= 100" --filter-name "Qual100" \
//     -filter "QD < 2.0" --filter-name "QD2" \
//     -filter "MQ < 35.0" --filter-name "MQ35" \
//     -filter "FS > 60.0" --filter-name "FS60" \
//     -filter "HaplotypeScore > 13.0" --filter-name "Haplo13" \
//     -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12" \
//     -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRank-8" \
//     -V ${vcf} -O ${scaffold}.markfiltered.vcf

//     gatk SelectVariants --exclude-filtered \
//       -V ${scaffold}.markfiltered.vcf \
//       -O ${scaffold}.filtered.vcf
// 	"""
// }

// // NB GatherVcfs must be provided with the vcf files in order of the scaffolds. We do this using the scaffhold_list.
// process gather_vcfs{
// 	tag "GatherVcfs"
//     container 'broadinstitute/gatk:4.1.9.0'

//     input:
//     val(scaffhold_list) from scaffold_list_gather_vcfs_ch.collect()
//     path(vcf) from vcf_gather_vcfs_ch.collect()
// 	path(vcf_idx) from vcf_idx_gather_vcfs_ch.collect()

// 	output:
//     tuple path("fowl.filtered.vcf"), path("fowl.filtered.vcf.idx") into make_bqsr_tables_known_variants_ch,make_bqsr_tables_known_variants_round_2_ch

//     script:
// 	"""
// 	gatk GatherVcfs ${scaffhold_list.collect{ "--INPUT ${it}.filtered.vcf " }.join()} --OUTPUT fowl.filtered.vcf
// 	"""
// }


// // We will want to do this on a per sample basis
// // We will want to use the fowl.vcf and the fowl.vcf.idx multiple times withhout using them up
// // It may be that we have to output these files as values rather than path/files in the gather_vcf process.
// // make_bqsr_tables_known_variants_ch.view()
// process make_bqsr_tables{
//     tag {pair_id}
//     container 'broadinstitute/gatk:4.1.9.0'

//     input:
//     tuple val(pair_id), path(merged_bam), path(known_variants), path(known_variants_idx) from make_bqsr_tables_bam_ch.combine(make_bqsr_tables_known_variants_ch)

//     output:
//     tuple val(pair_id), path(merged_bam), path("${pair_id}.table") into apply_bqsr_tables_ch
//     tuple val(pair_id), path("${pair_id}.table") into compare_tables_round_1_ch

//     script:
//     """
//     gatk BaseRecalibrator -I ${merged_bam[0]} -R ${params.ref_assembly_path} \
//     --known-sites $known_variants -O ${pair_id}.table
//     """
// }

// process apply_bqsr_tables{
//     tag {pair_id}
//     container 'broadinstitute/gatk:4.1.9.0'

//     input:
//     tuple val(pair_id), path(merged_bam), path(table) from apply_bqsr_tables_ch

//     output:
//     tuple val(pair_id), path("${pair_id}.recalibrated.bam{,.bai}") into apply_bqsr_tables_out_ch,apply_bqsr_tables_round_2_ch

//     script:
//     """
//     gatk ApplyBQSR -R ${params.ref_assembly_path} -I ${merged_bam[0]} \
//     --bqsr-recal-file $table \
//     -O ${pair_id}.recalibrated.bam;
//     """
// }

// // See the below issue to get an idea of the general strategy of getting AnalyzeCovariates to work in 
// // GATK 4.
// // https://github.com/broadinstitute/gatk/issues/322

// process make_bqsr_tables_round_2{
//     tag {pair_id}
//     container 'broadinstitute/gatk:4.1.9.0'

//     input:
//     tuple val(pair_id), path(merged_bam), path(known_variants), path(known_variants_idx) from apply_bqsr_tables_round_2_ch.combine(make_bqsr_tables_known_variants_round_2_ch)

//     output:
//     tuple val(pair_id), path("${pair_id}.after.table") into compare_tables_round_2_ch

//     script:
//     """
//     gatk BaseRecalibrator -I ${merged_bam[0]} -R ${params.ref_assembly_path} \
//     --known-sites $known_variants -O ${pair_id}.after.table
//     """

// }

// process AnalyzeCovariates{
//     tag {pair_id}
//     container 'broadinstitute/gatk:4.1.9.0'
//     publishDir analyze_covariates_publishDir

//     input:
//     tuple val(pair_id), path(table_before), path(table_after) from compare_tables_round_1_ch.join(compare_tables_round_2_ch)

//     output:
//     tuple val(pair_id), path("${pair_id}.BQSR.csv"), path("${pair_id}.BQSR.pdf") into compare_tables_round_2_out_ch
    
//     script:
//     """
//     gatk AnalyzeCovariates \
//       -before $table_before \
//       -after $table_after \
//       -csv ${pair_id}.BQSR.csv \
//       -plots ${pair_id}.BQSR.pdf
//     """
// }