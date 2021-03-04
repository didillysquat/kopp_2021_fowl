#!/usr/bin/env nextflow

// Check that the ref genome is decompressed
if (params.ref_assembly_path.endsWith(".bgz") || params.ref_assembly_path.endsWith(".gz") || params.ref_assembly_path.endsWith(".zip")){
    throw new Exception("The reference assembly genome must be decompressed and in fasta format (.fna, .fa, .fasta). Currently the reference genome path is set to ${params.ref_assembly_path}")
}

bin_dir = "${workflow.launchDir}/bin"
envs_dir = "${workflow.launchDir}/envs"
launch_dir = "${workflow.launchDir}"
tru_seq_pe_fasta_path = "${workflow.launchDir}/TruSeq3-PE.fa"


scaffold_list = {  
    def scaffolds = []
    new File(params.ref_assembly_path).eachLine {
        line -> 
        if (line.startsWith(">")){scaffolds << line.split()[0][1..-1];}
        }
    return scaffolds
}()

Channel.fromList(scaffold_list).into{scaffold_list_ch; scaffold_list_gather_vcfs_ch}

// This creates the read_group_map that is a list of tuples where the first element is the pair_id
// And the second element is the string that AddOrReplaceReadGroups will use to add readgroups and samples
// A unique sample name is given to each pair_id indepent of read group.
// Each of the files listed in the params.path_to_read_group_lists will get the same readgroup.
// This read_group_map is joined with the read_groups_merged_ch.
// For this to work the read pairs must be in the form <pair_id>_R[1,2].*.gz.
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
    params.output_dir = "${workflow.launchDir}/sub_sampled_outputs/pre_processing"
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim_sub_${params.subsample_depth}"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim_sub_${params.subsample_depth}"].join(File.separator)
    markduplicates_metrics_publishDir = [params.output_dir, "markduplicates_metrics_sub_${params.subsample_depth}"].join(File.separator)
    output_bam_publishDir = [params.output_dir, "output_bams_sub_${params.subsample_depth}"].join(File.separator)
    pre_seq_c_curve_publish_dir = [params.output_dir, "pre_seq_c_curve_sub_${params.subsample_depth}"].join(File.separator)
    collect_gc_bias_metrics_publishDir = [params.output_dir, "collect_gc_bias_metrics_sub_${params.subsample_depth}"].join(File.separator)
    pcr_bottleneck_coefficient_publishDir = [params.output_dir, "pcr_bottleneck_coefficient_sub_${params.subsample_depth}"].join(File.separator)
    mosdepth_sequencing_coverage_publishDir = [params.output_dir, "mosdepth_sequencing_coverage_sub_${params.subsample_depth}"].join(File.separator)
}else{
    params.output_dir = "${workflow.launchDir}/outputs/pre_processing"
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim"].join(File.separator)
    markduplicates_metrics_publishDir = [params.output_dir, "markduplicates_metrics"].join(File.separator)
    output_bam_publishDir = [params.output_dir, "output_bams"].join(File.separator)
    pre_seq_c_curve_publish_dir = [params.output_dir, "pre_seq_c_curve"].join(File.separator)
    collect_gc_bias_metrics_publishDir = [params.output_dir, "collect_gc_bias_metrics"].join(File.separator)
    pcr_bottleneck_coefficient_publishDir = [params.output_dir, "pcr_bottleneck_coefficient"].join(File.separator)
    mosdepth_sequencing_coverage_publishDir = [params.output_dir, "mosdepth_sequencing_coverage"].join(File.separator)
}

// START OF PRE PROCESSING
/* 
If subsample, create a channel that will pass into the subsampling process and then pass the subsampled
files into the trimming and fastqc
Modify the pair name to indicate subsampled files. This way it will carry through the remainder of the anlysis.
if not subsample, then create channel directly from the raw sequencing files
*/
if (params.subsample){
    ch_subsample = Channel.fromFilePairs("${params.raw_reads_dir}/*R{1,2}*.fastq.gz")    
    process subsample{
        tag "${pair_id}"
        cache 'lenient'
        container 'biocontainers/seqtk:v1.3-1-deb_cv1'
        containerOptions '-u $(id -u):$(id -g)'
        publishDir "${params.output_dir}/${params.subsample_depth}_subsampled_reads", mode: 'copy'
        cpus 1

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
}

process fastqc_pre_trim{
    tag "${pair_id}"
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir fastqc_pre_trim_publish_dir, mode: 'copy'
    cpus 1

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
    cpus params.trimmomatic_threads
    memory '24 GB'
	
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
    cpus 1
    memory '24 GB'

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
// Two options here, either supply the params.ref_assembly_path directly to ngm
// or path in ref_genome path made from the params.ref_assembly_path
// They have different results. For the first, ngm will look for the indexed files
// in the same directory as the params.ref_assembly_path. If they exist the indexing will be skipped.
// With this option the param arguments must be again passed directly to any process that uses the indices.
// Alternatively, with the second method, the indices will be stored in the nextflow work directory
// and we need to channel these appropriately into any process that wants to use them.
// The former relies on the params.ref_assembly_path being accessibly by next flow (i.e. below the root directory)
// of the project. The latter method does not. As such we will go with the latter.
process nextgenmap_indexing{
    conda 'envs/ngm.yaml'
    cpus params.nextgenmap_threads
    memory '24 GB'

    input:
    path ref_genome from params.ref_assembly_path

    output:
    tuple path(ref_genome), path("*.ngm") into ngm_index_ch

    script:
    """
    ngm -t ${task.cpus} -r ${ref_genome}
    """
}

process nextgenmap_mapping{
    tag "${pair_id}"
    conda 'envs/ngm.yaml'
    cpus params.nextgenmap_threads
    memory '24 GB'

    input:
	tuple val(pair_id), file(paired), file(unpaired), path(ref_genome), path(ref_genome_indices) from ch_ngm_paired.join(ch_ngm_unpaired).combine(ngm_index_ch)

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
    cpus 1
    memory '24 GB'

    input:
    tuple val(pair_id), file(paired), file(unpaired_one), file(unpaired_two) from add_read_group_headers_ch

    output:
    tuple val(pair_id), file("${pair_id}.merged.mapped.sam") into read_groups_merged_ch

    script:
    """
    gatk MergeSamFiles --INPUT $paired --INPUT $unpaired_one --INPUT $unpaired_one --OUTPUT ${pair_id}.merged.mapped.sam
    """

}

process add_read_group_headers{
    tag "${pair_id}"
    container 'broadinstitute/gatk:4.1.9.0'
    cpus 1
    memory '24 GB'

    input:
    tuple val(pair_id), file(merged), val(read_group_string) from read_groups_merged_ch.join(Channel.fromList(read_group_map))

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
    publishDir markduplicates_metrics_publishDir, pattern: "*.metrics.txt", mode: 'move'
    publishDir output_bam_publishDir, pattern: "*.bam{,.bai}", mode: 'copy'
    cpus 4
    memory "24 GB"

    input:
    tuple val(pair_id), file(merged) from mark_duplicates_ch

    output:
    tuple val(pair_id), file("${pair_id}.merged.deduplicated.sorted.bam{,.bai}") into collect_gc_bias_metrics_ch,pcr_bottleneck_coefficient_ch,mosdepth_sequencing_coverage_ch
    file("${pair_id}.merged.deduplicated.sorted.metrics.txt") into mark_duplicate_metrics_ch

    script:
    """
    gatk MarkDuplicatesSpark --create-output-bam-index --remove-sequencing-duplicates -I ${merged} -O ${pair_id}.merged.deduplicated.sorted.bam -M ${pair_id}.merged.deduplicated.sorted.metrics.txt
    """
}
// END OF PRE PROCESSING

// // START OF EVALUATION METRICS
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