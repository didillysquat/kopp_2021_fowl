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

make_indices_from_scratch = {
    if (!params.remake_indices){
        // Check to see if the reference assembly indices already exist
        if (params.mapping == "ngm"){
            from_scratch = {   
                def file_objects = [];
                def extensions = ["", "-enc.2.ngm", "-ht-13-2.3.ngm"];
                def all_indices_found = true;
                extensions.each{
                    def file_object = new File("${params.ref_assembly_path}${it}")
                    if (!file_object.exists()) {all_indices_found = false;}
                }
                return !all_indices_found
            }()
            return from_scratch
        }else{
            from_scratch = {   
                def file_objects = [];
                def extensions = ["", ".sa", ".amb", ".ann", ".pac", ".bwt"];
                def all_indices_found = true;
                extensions.each{
                    def file_object = new File("${params.ref_assembly_path}${it}")
                    if (!file_object.exists()) {all_indices_found = false;}
                }
                return !all_indices_found
            }()
            return from_scratch
        }
    }else{
        // Recompute the reference assmbly indices irrespective of if they already exist
        return true
    }
}()



// This method makes read_group_map that is a list of tuples where the first element is the pair_id
// and the second is the string that is used by AddOrReplaceReadGroups.

(read_group_map, num_samples) = {
    def tup_list = []
    def pair_id_list = []
    new File(params.read_info_file).eachLine {  
        line ->
        line_comp = line.tokenize();
        // Skip the header line of the .tsv
        if (line_comp[0].endsWith(".gz")){            
            // TODO check to see that the file exists and return error if it does not
            def should_exist = new File("${params.raw_reads_dir}/${line_comp[0]}")
            if (!should_exist.exists()){
                throw new Exception("${should_exist} not found.")
            }
            
            // regex to get the pair_id
            def pattern = ~/^(?<pair>.*)_R[1-2].*$/
            def matcher = line_comp[0] =~ pattern
            
            if (matcher){
                pair_id = matcher.group("pair")
                if (params.subsample){
                    pair_id = "${pair_id}_sub_${params.subsample_depth}"
                }
                // At this point we need to see if the mapping information has already been added
                // If not, then add in the mapping information 
                if (!pair_id_list.contains(pair_id)){
                    pair_id_list << pair_id
                    tup_list << ["${pair_id}".toString(), "RGID=${line_comp[1]} RGLB=${line_comp[2]} RGPL=${line_comp[3]} RGPU=${line_comp[4]} RGSM=${line_comp[5]}"]
                }                            
            }else{
                throw new Exception("An error has occured when making the readgroup dictionary\nNo match was found for ${line_comp[0]}")
            }
        }
    }
    return [tup_list, pair_id_list.size]
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
// NB the containerOoptions are required when using the docker profile to prevent permission errors
// They are not required when running a singularity profile.
if (params.subsample){
    ch_subsample = Channel.fromFilePairs("${params.raw_reads_dir}/*R{1,2}*.fastq.gz")    
    process subsample{
        tag "${pair_id}"
        cache 'lenient'
        container 'biocontainers/seqtk:v1.3-1-deb_cv1'
        if (workflow.containerEngine == 'docker'){
            containerOptions '-u $(id -u):$(id -g)'
        }
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
    Channel.fromFilePairs("${params.raw_reads_dir}/*R{1,2}*.fastq.gz").into{ch_fastqc_pre_trim; ch_trimmomatic_input}
}

// TODO these first three processes should likely be replaced with fastp in future pipelines.
// NB the biocontainers container biocontainers/fastqc:v0.11.9_cv7 does not work on BINAC.
process fastqc_pre_trim{
    tag "${pair_id}"
    container 'singlecellpipeline/fastqc:v0.0.2'
    publishDir fastqc_pre_trim_publish_dir, mode: 'copy'
    cpus 1

    input:
    tuple val(pair_id), file(fastq_file) from ch_fastqc_pre_trim.flatMap{[["${it[0]}_1", it[1][0]], ["${it[0]}_2", it[1][1]]]}

    output:
    file "${out_name}" into ch_fastqc_pre_trim_output

    script:
    out_name = fastq_file.getName().replaceAll('.fastq.gz', '.pre_trim.fastqc.html')
    """
    fastqc -t 2 -o . $fastq_file
    mv *.html ${out_name}
    """
}

process trimmomatic{
    cache 'lenient'
	tag "${pair_id}"
    container 'davelabhub/trimmomatic:0.39--1'
    if (workflow.containerEngine == 'docker'){
            containerOptions '-u $(id -u):$(id -g)'
        }
    cpus params.trimmomatic_threads
	
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

// NB the biocontainers container biocontainers/fastqc:v0.11.9_cv7 does not work on BINAC.
// NB running fastqc on the unpaired files is causing fastqc to hang so we will only run
// on the paired files: https://github.com/s-andrews/FastQC/issues/74
// As an extra causion we will use the -t 2 option to use two threads wich each get allocated
// 250 MB of memory.
process fastqc_post_trim{
    tag "${pair_id}"
    container 'singlecellpipeline/fastqc:v0.0.2'
    publishDir fastqc_post_trim_publish_dir, mode: 'copy'
    cpus 1

    input:
    tuple val(pair_id), file(paired_files), file(unpaired_files) from ch_fastqc_post_trim_paired.join(ch_fastqc_post_trim_unpaired)

    output:
    tuple file("${pair_id}*1P*fastqc.html"), file("${pair_id}*2P*fastqc.html") into ch_fastqc_post_trim_output

    script:
    """
    fastqc -t 2 -o . ${paired_files[0]}
    fastqc -t 2 -o . ${paired_files[1]}
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

// NB we check to see if the indexing files already exist and
// make use of them if params.remake_indices if false
if (params.mapping == "ngm"){
    if (make_indices_from_scratch){
        process nextgenmap_indexing{
            conda 'envs/ngm.yaml'
            cpus params.nextgenmap_threads

            input:
            path ref_genome from params.ref_assembly_path

            output:
            tuple path(ref_genome), path("*.ngm") into ngm_index_ch

            script:
            """
            ngm -t ${task.cpus} -r ${ref_genome}
            """
        }
    }else{
        // Manually create the channel
        println("Using existing NextGenMap indices")
        ngm_index_ch = Channel.fromList([[file(params.ref_assembly_path), file("${params.ref_assembly_path}{-enc.2.ngm, -ht-13-2.3.ngm}")]])
    }
    
    process nextgenmap_mapping{
        tag "${pair_id}"
        conda 'envs/ngm.yaml'
        cpus params.nextgenmap_threads

        input:
        tuple val(pair_id), file(paired), file(unpaired), path(ref_genome), path(ref_genome_indices) from ch_ngm_paired.join(ch_ngm_unpaired).combine(ngm_index_ch)

        output:
        tuple val(pair_id), file("${pair_id}_P_mapped.sam"), file("${pair_id}_U1_mapped.sam"), file("${pair_id}_U2_mapped.sam") into add_read_group_headers_ch

        script:
        """
        ngm -t ${task.cpus} -r ${ref_genome} -1 ${paired[0]} -2 ${paired[1]} -o ${pair_id}_P_mapped.sam;
        ngm -t ${task.cpus} -r ${ref_genome} -q ${unpaired[0]} -o ${pair_id}_U1_mapped.sam;
        ngm -t ${task.cpus} -r ${ref_genome} -q ${unpaired[1]} -o ${pair_id}_U2_mapped.sam;
        """
    }
}else{
    // Check to see if the indice files already exist and if create the channel from this
    if (make_indices_from_scratch){
        process bwa_indexing{
            container 'biocontainers/bwa:v0.7.17_cv1'
            if (workflow.containerEngine == 'docker'){
                containerOptions '-u $(id -u):$(id -g)'
            }
            cpus params.nextgenmap_threads

            input:
            path ref_genome from params.ref_assembly_path

            output:
            tuple path(ref_genome), path("*{.sa,.amb,.ann,.pac,.bwt}") into ngm_index_ch

            script:
            """
            bwa index ${ref_genome}
            """
        }
    }else{
        // Manually create the channel 
        println("Using existing BWA indices")
        ngm_index_ch = Channel.fromList([[file(params.ref_assembly_path), file("${params.ref_assembly_path}{.sa,.amb,.ann,.pac,.bwt}")]])
    }

    process bwa_mapping{
        tag "${pair_id}"
        container 'biocontainers/bwa:v0.7.17_cv1'
        if (workflow.containerEngine == 'docker'){
            containerOptions '-u $(id -u):$(id -g)'
        }
        cpus params.nextgenmap_threads

        input:
        tuple val(pair_id), file(paired), file(unpaired), path(ref_genome), path(ref_genome_indices) from ch_ngm_paired.join(ch_ngm_unpaired).combine(ngm_index_ch)

        output:
        tuple val(pair_id), file("${pair_id}_P_mapped.unsorted.sam"), file("${pair_id}_U1_mapped.unsorted.sam"), file("${pair_id}_U2_mapped.unsorted.sam") into samtools_sort_ch

        script:
        """
        bwa mem -t ${task.cpus} ${ref_genome} ${paired[0]} ${paired[1]} > ${pair_id}_P_mapped.unsorted.sam;
        bwa mem -t ${task.cpus} ${ref_genome} ${unpaired[0]} > ${pair_id}_U1_mapped.unsorted.sam;
        bwa mem -t ${task.cpus} ${ref_genome} ${unpaired[1]} > ${pair_id}_U2_mapped.unsorted.sam;
        """
    }
    
    // --threads is the number of additional threads to use (hence the -1)
    process samtools_sort{
        tag "${pair_id}"
        container 'singlecellpipeline/samtools:v0.0.3'
        cpus params.nextgenmap_threads

        input:
        tuple val(pair_id), file(paired), file(unpaired_1), path(unpaired_2) from samtools_sort_ch

        output:
        tuple val(pair_id), file("${pair_id}_P_mapped.sam"), file("${pair_id}_U1_mapped.sam"), file("${pair_id}_U2_mapped.sam") into add_read_group_headers_ch

        script:
        """
        samtools sort --threads ${task.cpus - 1} ${paired} --output-fmt SAM > ${pair_id}_P_mapped.sam;
        samtools sort --threads ${task.cpus - 1} ${unpaired_1} --output-fmt SAM > ${pair_id}_U1_mapped.sam;
        samtools sort --threads ${task.cpus - 1} ${unpaired_2} --output-fmt SAM > ${pair_id}_U2_mapped.sam;
        """
    }
}

process merge_paired_and_unpaired{
    tag "${pair_id}"
    container 'broadinstitute/gatk:4.2.0.0'
    cpus 1

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
    container 'broadinstitute/gatk:4.2.0.0'
    cpus 1

    input:
    tuple val(pair_id), file(merged), val(read_group_string) from read_groups_merged_ch.join(Channel.fromList(read_group_map))

    output:
    tuple val(pair_id), file("${pair_id}.merged.mapped.headers.sam") into mark_duplicates_ch

    script:
    """
    gatk AddOrReplaceReadGroups I=${merged} O=${pair_id}.merged.mapped.headers.sam ${read_group_string}
    """
}

// NB when running a singularity profile we get errors relating to port number conflicts.
// To prevent this we set spark.port.maxRetries to be the number of samples
// This does not appear to be an issue when running a docker profile.
process markduplicates_spark{
    tag "${pair_id}"
    container 'broadinstitute/gatk:4.2.0.0'
    publishDir markduplicates_metrics_publishDir, pattern: "*.metrics.txt", mode: 'copy'
    publishDir output_bam_publishDir, pattern: "*.bam{,.bai}", mode: 'copy'
    cpus 4

    input:
    tuple val(pair_id), file(merged) from mark_duplicates_ch

    output:
    tuple val(pair_id), file("${pair_id}.merged.deduplicated.sorted.bam{,.bai}") into collect_gc_bias_metrics_ch,pcr_bottleneck_coefficient_ch,mosdepth_sequencing_coverage_ch
    file("${pair_id}.merged.deduplicated.sorted.metrics.txt") into mark_duplicate_metrics_ch

    script:
    """
    gatk MarkDuplicatesSpark --create-output-bam-index --remove-sequencing-duplicates -I ${merged} -O ${pair_id}.merged.deduplicated.sorted.bam -M ${pair_id}.merged.deduplicated.sorted.metrics.txt --conf 'spark.port.maxRetries=${num_samples}'
    """
}
// END OF PRE PROCESSING

// // START OF EVALUATION METRICS
process collect_gc_bias_metrics{
    tag "${pair_id}"
    publishDir collect_gc_bias_metrics_publishDir, mode: 'copy'
    container 'broadinstitute/gatk:4.2.0.0'
    cpus 1

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
    tag "${pair_id}"
    publishDir pcr_bottleneck_coefficient_publishDir, mode: 'copy'
    container 'encodedcc/atac-seq-pipeline:PIP-1469_pbam_b92239a0-82a0-4297-b8b3-3a655a4626a8'
    cpus 1

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
    tag "${pair_id}"
    publishDir mosdepth_sequencing_coverage_publishDir, mode: 'copy'
    container 'davelabhub/mosdepth:0.2.5--hb763d49_0'
    cpus 1
    
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
    tag "${pair_id}"
    publishDir mosdepth_sequencing_coverage_publishDir, mode: 'copy'
    container "quay.io/biocontainers/python:3.8.3"
    cpus 1

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