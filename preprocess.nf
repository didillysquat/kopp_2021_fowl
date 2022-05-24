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

/* TODO
Make the input format come from tsv specifiying full sequence paths, read group info and sample name
Adapt to work with SE and PE reads
Add the pileupCaller
Add FreeBayes
Check that subsampling works by only providing the --subsample flag
*/

// Check that the ref genome is decompressed
if (params.ref.endsWith(".bgz") || params.ref.endsWith(".gz") || params.ref.endsWith(".zip")){
    throw new Exception("The reference assembly genome must be decompressed and in fasta format (.fna, .fa, .fasta). Currently the reference genome path is set to ${params.ref}")
}

bin_dir = "${workflow.launchDir}/bin"
envs_dir = "${workflow.launchDir}/envs"
launch_dir = "${workflow.launchDir}"
tru_seq_pe_fasta_path = "${workflow.launchDir}/TruSeq3-PE.fa"
def tru_seq_file = new File("${tru_seq_pe_fasta_path}")
if (!tru_seq_file.exists()) {throw new Exception("Could not find TruSeq3-PE.fa")}

if (!params.trimmomatic_threads){
    params.trimmomatic_threads = 1
}
if (!params.mapping_threads){
    params.mapping_threads = 1
}
if (!params.mark_duplicate_threads){
    params.mark_duplicate_threads = 1
}

if (!params.estlibcomp_max_mem){
    params.estlibcomp_max_mem = 8
}

params.remake_indices = false;

if (!params.subsample){
    params.subsample = false;
}

scaffold_list = {  
    def scaffolds = []
    new File(params.ref).eachLine {
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
                    def file_object = new File("${params.ref}${it}")
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
                    def file_object = new File("${params.ref}${it}")
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
// It also outputs the number of samples and sets up the initial Channel that holds the reads

(read_group_map, num_samples, reads_ch) = {
    def tup_list = []
    def pair_id_list = []
    def reads_ch_list = []
    new File(params.input_tsv).eachLine {  
        line ->
        line_comp = line.tokenize();
        // Skip the header line of the .tsv
        if (line_comp[0].endsWith(".gz")){            
            // Check to see that the file exists and return error if it does not
            def should_exist_one = new File(line_comp[0])
            def should_exist_two = new File(line_comp[1])
            if (!should_exist_one.exists()){
                throw new Exception("${should_exist_one} not found.")
            }
            if (!should_exist_two.exists()){
                throw new Exception("${should_exist_two} not found.")
            }
            
            pair_id = line_comp[6]
            if (params.subsample){
                    pair_id = "${pair_id}_sub_${params.subsample_depth}"
                }
            // At this point we need to see if the mapping information has already been added
            // If not, then add in the mapping information 
            if (!pair_id_list.contains(pair_id)){
                pair_id_list << pair_id
                tup_list << ["${pair_id}".toString(), "--RGID ${line_comp[2]} --RGLB ${line_comp[3]} --RGPL ${line_comp[4]} --RGPU ${line_comp[5]} --RGSM ${line_comp[6]}"]
            }else{
                throw new Exception("There appear to be non-uniue sample names in your input.tsv: ${pair_id}")
            }
            reads_ch_list << [pair_id, [line_comp[0], line_comp[1]]]
        }
    }
    return [tup_list, pair_id_list.size, reads_ch_list]
}()

// Publish dirs
if (params.subsample){
    params.output_dir = "${workflow.launchDir}/sub_sampled_outputs/pre_processing"
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim_sub_${params.subsample_depth}"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim_sub_${params.subsample_depth}"].join(File.separator)
    markduplicates_metrics_publishDir = [params.output_dir, "markduplicates_metrics_sub_${params.subsample_depth}"].join(File.separator)
    estimatelibrarycomplexity_metrics_publishDir = [params.output_dir, "estimatelibrarycomplexity_metrics_sub_${params.subsample_depth}"].join(File.separator)
    output_bam_publishDir = [params.output_dir, "output_bams_sub_${params.subsample_depth}"].join(File.separator)
    preseq_complexity_prediction_publishDir = [params.output_dir, "preseq_complexity_sub_${params.subsample_depth}"].join(File.separator)
    collect_gc_bias_metrics_publishDir = [params.output_dir, "collect_gc_bias_metrics_sub_${params.subsample_depth}"].join(File.separator)
    pcr_bottleneck_coefficient_publishDir = [params.output_dir, "pcr_bottleneck_coefficient_sub_${params.subsample_depth}"].join(File.separator)
    mosdepth_sequencing_coverage_publishDir = [params.output_dir, "mosdepth_sequencing_coverage_sub_${params.subsample_depth}"].join(File.separator)
    samtools_coverage_stats_publishDir = [params.output_dir, "samtools_coverage_stats_sub_${params.subsample_depth}"].join(File.separator)
    samtools_mapping_stats_prededup_publishDir = [params.output_dir, "samtools_mapping_stats_prededup_sub_${params.subsample_depth}"].join(File.separator)
    samtools_mapping_stats_postdedup_publishDir = [params.output_dir, "samtools_mapping_stats_postdedup_sub_${params.subsample_depth}"].join(File.separator)
    preprocessing_summary_publishDir = [params.output_dir, "preprocessing_summary_sub_${params.subsample_depth}"].join(File.separator)
}else{
    params.output_dir = "${workflow.launchDir}/outputs/pre_processing"
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim"].join(File.separator)
    markduplicates_metrics_publishDir = [params.output_dir, "markduplicates_metrics"].join(File.separator)
    estimatelibrarycomplexity_metrics_publishDir = [params.output_dir, "estimatelibrarycomplexity_metrics"].join(File.separator)
    output_bam_publishDir = [params.output_dir, "output_bams"].join(File.separator)
    preseq_complexity_prediction_publishDir = [params.output_dir, "preseq_complexity"].join(File.separator)
    collect_gc_bias_metrics_publishDir = [params.output_dir, "collect_gc_bias_metrics"].join(File.separator)
    pcr_bottleneck_coefficient_publishDir = [params.output_dir, "pcr_bottleneck_coefficient"].join(File.separator)
    mosdepth_sequencing_coverage_publishDir = [params.output_dir, "mosdepth_sequencing_coverage"].join(File.separator)
    samtools_coverage_stats_publishDir = [params.output_dir, "samtools_coverage_stats"].join(File.separator)
    samtools_mapping_stats_prededup_publishDir = [params.output_dir, "samtools_mapping_stats_prededup"].join(File.separator)
    samtools_mapping_stats_postdedup_publishDir = [params.output_dir, "samtools_mapping_stats_postdedup"].join(File.separator)
    preprocessing_summary_publishDir = [params.output_dir, "preprocessing_summary"].join(File.separator)
}

input_tsv = file(params.input_tsv)
// START OF PRE PROCESSING
/* 
If subsample, create a channel that will pass into the subsampling process and then pass the subsampled
files into the trimming and fastqc
Modify the pair name to indicate subsampled files. This way it will carry through the remainder of the anlysis.
if not subsample, then create channel directly from the raw sequencing files
*/
// NB the container options are required when using the docker profile to prevent permission errors
// They are not required when running a singularity profile.
if (params.subsample){
    ch_subsample = Channel.fromList(reads_ch)
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
        tuple val(pair_id), path(reads) from ch_subsample

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
    Channel.fromList(reads_ch).into{ch_fastqc_pre_trim; ch_trimmomatic_input}
}


// NB the biocontainers container biocontainers/fastqc:v0.11.9_cv7 does not work on BINAC.
process fastqc_pre_trim{
    tag "${pair_id}"
    container 'singlecellpipeline/fastqc:v0.0.2'
    publishDir fastqc_pre_trim_publish_dir, mode: 'copy'
    cpus 1
    cache 'lenient'

    input:
    tuple val(pair_id), path(fastq_file) from ch_fastqc_pre_trim.flatMap{[["${it[0]}_1", it[1][0]], ["${it[0]}_2", it[1][1]]]}

    output:
    file "${out_name}" into ch_fastqc_pre_trim_output

    script:
    out_name = fastq_file.getName().replaceAll('.fastq.gz', '.pre_trim.fastqc.html').replaceAll('.fq.gz', '.pre_trim.fastqc.html')
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
	tuple val(pair_id), path(fastqs) from ch_trimmomatic_input
	
	output:
	tuple val(pair_id), file("${pair_id}*{1,2}P.fq.gz") into ch_fastqc_post_trim_paired,ch_ngm_paired
    tuple val(pair_id), file("${pair_id}*{1,2}U.fq.gz") into ch_fastqc_post_trim_unpaired,ch_ngm_unpaired

	script:
	outbase = fastqs[0].getName().replaceAll('.fastq.gz', '.trimmed.fq.gz').replaceAll('.fq.gz', '.trimmed.fq.gz')
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
// As an extra caution we will use the -t 2 option to use two threads wich each get allocated
// 250 MB of memory.
process fastqc_post_trim{
    tag "${pair_id}"
    container 'singlecellpipeline/fastqc:v0.0.2'
    publishDir fastqc_post_trim_publish_dir, mode: 'copy'
    cpus 1
    cache 'lenient'

    input:
    tuple val(pair_id), file(paired_files), file(unpaired_files) from ch_fastqc_post_trim_paired.join(ch_fastqc_post_trim_unpaired)

    output:
    tuple file("${pair_id}*1P*fastqc.html"), file("${pair_id}*2P*fastqc.html"), file("${pair_id}*1U*fastqc.html"), file("${pair_id}*2U*fastqc.html") into ch_fastqc_post_trim_output

    script:
    """
    fastqc -t 2 -o . ${paired_files[0]}
    fastqc -t 2 -o . ${paired_files[1]}
    fastqc -t 2 -o . ${unpaired_files[0]}
    fastqc -t 2 -o . ${unpaired_files[1]}
    """
}

// Index ref assembly and map
// NB we check to see if the indexing files already exist and
// make use of them if params.remake_indices if false
if (params.mapping == "ngm"){
    if (make_indices_from_scratch){
        process nextgenmap_indexing{
            conda 'envs/ngm.yaml'
            cpus params.mapping_threads

            input:
            path ref_genome from params.ref

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
        ngm_index_ch = Channel.fromList([[file(params.ref), file("${params.ref}{-enc.2.ngm, -ht-13-2.3.ngm}")]])
    }
    
    process nextgenmap_mapping{
        tag "${pair_id}"
        conda 'envs/ngm.yaml'
        cpus params.mapping_threads

        input:
        tuple val(pair_id), file(paired), file(unpaired), path(ref_genome), path(ref_genome_indices) from ch_ngm_paired.join(ch_ngm_unpaired).combine(ngm_index_ch)

        output:
        tuple val(pair_id), file("${pair_id}_P_mapped.unsorted.bam"), file("${pair_id}_U1_mapped.unsorted.bam"), file("${pair_id}_U2_mapped.unsorted.bam") into samtools_sort_ch

        script:
        """
        ngm -b -t ${task.cpus} -r ${ref_genome} -1 ${paired[0]} -2 ${paired[1]} -o ${pair_id}_P_mapped.unsorted.bam;
        ngm -b -t ${task.cpus} -r ${ref_genome} -q ${unpaired[0]} -o ${pair_id}_U1_mapped.unsorted.bam;
        ngm -b -t ${task.cpus} -r ${ref_genome} -q ${unpaired[1]} -o ${pair_id}_U2_mapped.unsorted.bam;
        """
    }
   
    process samtools_sort{
        tag "${pair_id}"
        container 'singlecellpipeline/samtools:v0.0.3'
        cpus params.mapping_threads

        input:
        tuple val(pair_id), file(paired), file(unpaired_1), path(unpaired_2) from samtools_sort_ch

        output:
        tuple val(pair_id), file("${pair_id}_P_mapped.bam"), file("${pair_id}_U1_mapped.bam"), file("${pair_id}_U2_mapped.bam") into merge_paired_unpaired_ch

        script:
        """
        samtools sort --threads ${task.cpus - 1} ${paired} --output-fmt BAM > ${pair_id}_P_mapped.bam;
        samtools sort --threads ${task.cpus - 1} ${unpaired_1} --output-fmt BAM > ${pair_id}_U1_mapped.bam;
        samtools sort --threads ${task.cpus - 1} ${unpaired_2} --output-fmt BAM > ${pair_id}_U2_mapped.bam;
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
            cpus params.mapping_threads

            input:
            path ref_genome from params.ref

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
        ngm_index_ch = Channel.fromList([[file(params.ref), file("${params.ref}{.sa,.amb,.ann,.pac,.bwt}")]])
    }

    process bwa_mapping_sort{
        tag "${pair_id}"
        container 'dukegcb/bwa-samtools:latest'
        if (workflow.containerEngine == 'docker'){
            containerOptions '-u $(id -u):$(id -g)'
        }
        cpus params.mapping_threads

        input:
        tuple val(pair_id), file(paired), file(unpaired), path(ref_genome), path(ref_genome_indices) from ch_ngm_paired.join(ch_ngm_unpaired).combine(ngm_index_ch)

        output:
        tuple val(pair_id), file("${pair_id}_P_mapped.bam"), file("${pair_id}_U1_mapped.bam"), file("${pair_id}_U2_mapped.bam") into merge_paired_unpaired_ch

        script:
        """
        bwa mem -t ${task.cpus} ${ref_genome} ${paired[0]} ${paired[1]} | samtools sort -@${task.cpus} -o ${pair_id}_P_mapped.bam -
        bwa mem -t ${task.cpus} ${ref_genome} ${unpaired[0]} | samtools sort -@${task.cpus} -o ${pair_id}_U1_mapped.bam -
        bwa mem -t ${task.cpus} ${ref_genome} ${unpaired[1]} | samtools sort -@${task.cpus} -o ${pair_id}_U2_mapped.bam -
        """
    }
    
}

process merge_paired_and_unpaired{
    tag "${pair_id}"
    container 'broadinstitute/gatk:4.2.4.1'
    cpus 1

    input:
    tuple val(pair_id), file(paired), file(unpaired_one), file(unpaired_two) from merge_paired_unpaired_ch

    output:
    tuple val(pair_id), file("${pair_id}.merged.mapped.bam"), file("${pair_id}.merged.mapped*.bai") into add_read_group_headers_ch

    script:
    """
    gatk --java-options "-XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" MergeSamFiles --SORT_ORDER coordinate --CREATE_INDEX --INPUT $paired --INPUT $unpaired_one --INPUT $unpaired_two --OUTPUT ${pair_id}.merged.mapped.bam
    """

}

process add_read_group_headers{
    tag "${pair_id}"
    container 'broadinstitute/gatk:4.2.4.1'
    cpus 1

    input:
    tuple val(pair_id), file(merged), file(merged_bai), val(read_group_string) from add_read_group_headers_ch.join(Channel.fromList(read_group_map))

    output:
    tuple val(pair_id), file("${pair_id}.merged.mapped.readGroupHeaders.bam"), file("${pair_id}.merged.mapped.readGroupHeaders*.bai") into mark_duplicates_ch, mapping_stats_prededup_ch, estimate_library_complexity_ch, bamtobed_ch

    script:
    """
    gatk --java-options "-XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" AddOrReplaceReadGroups --CREATE_INDEX true --INPUT ${merged} --OUTPUT ${pair_id}.merged.mapped.readGroupHeaders.bam ${read_group_string}
    """
}

// NB the preseq c_curve and lc_extrap functions seem to be broken for
// pe bam files. A work around is to convert the bam to a bed and supply
// the bed as input to the two preseq functions.
// I have made a container to do this in one process so that we can delete
// the intermediary bed file as these are about the same size as the original
// bam.
// 20220215 see https://github.com/smithlabcode/preseq/issues/59
// https://github.com/smithlabcode/preseq/issues/57
process bamtobed_preseq{
    tag "${pair_id}"
    container "didillysquat/bedtools_preseq:latest"
    publishDir preseq_complexity_prediction_publishDir, mode: "copy"

    input:
    tuple val(pair_id), file(bam), file(bai) from bamtobed_ch

    output:
    tuple file("${pair_id}.c_curve.txt"), file("${pair_id}.lc_extrap.txt") into preseq_out_ch

    script:
    """
    bedtools bamtobed -i $bam > ${pair_id}.merged.mapped.readGroupHeaders.bed
    preseq c_curve -P ${pair_id}.merged.mapped.readGroupHeaders.bed > ${pair_id}.c_curve.txt
    preseq lc_extrap -P ${pair_id}.merged.mapped.readGroupHeaders.bed > ${pair_id}.lc_extrap.txt
    rm ${pair_id}.merged.mapped.readGroupHeaders.bed
    """
}

// We will plot up a single figure that contains all of the 
// c_cures and lc_extrap cures in a single figure
// If users want a plot for each sample they will be able to make this from
// individual c_curve.txt and lc_extrap.txt files that were published from 
// the bamtobed_preseq process above.
process plot_preseq_complexity_curves{
    tag "plot_preseq_complexity_curves"
    container "didillysquat/r_multipurpose:latest"
    publishDir preseq_complexity_prediction_publishDir, mode: "copy"

    input:
    path curve_files from preseq_out_ch.collect()

    output:
    file "library.complexity.png" into plot_preseq_complexity_curves_out_ch

    script:
    """
    Rscript ${bin_dir}/plot_preseq_lib_complex.r \$PWD
    """

}


process samtools_mapping_stats_pre_deduplication{
    tag "${sample}"
    publishDir samtools_mapping_stats_prededup_publishDir, mode: "copy"
    container 'singlecellpipeline/samtools:v0.0.3'

    input:
    tuple val(sample), path(bam), path(bai) from mapping_stats_prededup_ch

    output:
    path "${sample}.bam.prededup.stats.txt" into samtools_stats_prededup_out

    script:
    """
    samtools stats $bam > ${sample}.bam.prededup.stats.txt
    """
}

// NB when running a singularity profile we get errors relating to port number conflicts.
// To prevent this we set spark.port.maxRetries to be the number of samples
// This does not appear to be an issue when running a docker profile.
// NB MarkDuplicatesSpark estimation of library complexity is based on only the mapped reads
// and is based on mapping position.
// By contrast EstimateLibraryComplexity is not based on mapping and so uses mapped and non-mapped
// paried reads. They can give very different results. We were in particular seeing
// that the MarkDuplicatesSpark estimate is strongly effected by the quality/proportion of reads
// mapping to reference genome. As such we will include both methods here and report them
// in the preprocessing summary output.
process markduplicates_spark{
    tag "${pair_id}"
    container 'broadinstitute/gatk:4.2.4.1'
    publishDir markduplicates_metrics_publishDir, pattern: "*.metrics.txt", mode: 'copy'
    publishDir output_bam_publishDir, pattern: "*.bam{,.bai}", mode: 'copy'
    cpus params.mark_duplicate_threads

    input:
    tuple val(pair_id), file(merged), file(merged_bai) from mark_duplicates_ch

    output:
    tuple val(pair_id), file("${pair_id}.merged.deduplicated.sorted.bam{,.bai}") into collect_gc_bias_metrics_ch,pcr_bottleneck_coefficient_ch,mosdepth_sequencing_coverage_ch,samtools_depth_stats_ch,samtools_mapping_stats_postdedup_ch,perprocess_overview_bams_ch
    file("${pair_id}.merged.deduplicated.sorted.metrics.txt") into mark_duplicate_metrics_ch

    script:
    """
    gatk MarkDuplicatesSpark --java-options "-XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" --create-output-bam-index true --remove-all-duplicates -I ${merged} -O ${pair_id}.merged.deduplicated.sorted.bam -M ${pair_id}.merged.deduplicated.sorted.metrics.txt --conf 'spark.port.maxRetries=${num_samples}' --spark-master local[${task.cpus}]
    """
}

process estimate_library_complexity{
    tag "${pair_id}"
    container 'broadinstitute/gatk:4.2.4.1'
    publishDir estimatelibrarycomplexity_metrics_publishDir, pattern: "*.metrics.txt", mode: 'copy'

    input:
    tuple val(pair_id), file(merged), file(merged_bai) from estimate_library_complexity_ch

    output:
    file("${pair_id}.estlibcomp.txt") into est_lib_comp_metrics_ch

    script:
    """
    gatk EstimateLibraryComplexity --java-options "-Xmx${params.estlibcomp_max_mem}g -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1" -I ${merged} -O ${pair_id}.estlibcomp.txt
    """
}

// END OF PRE PROCESSING

// // START OF EVALUATION METRICS
process collect_gc_bias_metrics{
    tag "${pair_id}"
    publishDir collect_gc_bias_metrics_publishDir, mode: 'copy'
    container 'broadinstitute/gatk:4.2.4.1'
    cpus 1

    input:
    tuple val(pair_id), path(merged) from collect_gc_bias_metrics_ch
    path ref_genome from params.ref
    
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
    maxRetries 5
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }

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

process samtools_mapping_stats_post_deduplication{
    tag "${sample}"
    publishDir samtools_mapping_stats_postdedup_publishDir, mode: "copy"
    container 'singlecellpipeline/samtools:v0.0.3'

    input:
    tuple val(sample), path(bam) from samtools_mapping_stats_postdedup_ch

    output:
    path "${sample}.bam.postdedup.stats.txt" into samtools_stats_postdedup_out

    script:
    """
    samtools stats ${bam[0]} > ${sample}.bam.postdedup.stats.txt
    """
}

process samtools_coverage_stats{
    tag "${sample}"
    publishDir samtools_coverage_stats_publishDir, mode: "copy"
    container 'singlecellpipeline/samtools:v0.0.3'

    input:
    tuple val(sample), path(bam) from samtools_depth_stats_ch

    output:
    path "${sample}.depth.txt" into samtools_depth_out

    script:
    """
    samtools depth -a ${bam[0]} | LC_NUMERIC=en_US.UTF-8 awk '{sum+=\$3; sumsq+=\$3*\$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)^2)}' > ${sample}.depth.txt
    """
}

process preprocess_summary{
    tag "preprocess_summary"
    publishDir preprocessing_summary_publishDir, mode: "copy"
    container 'amancevice/pandas:1.3.1'

    input:
    path pretrim from ch_fastqc_pre_trim_output.collect()
    path posttrim from ch_fastqc_post_trim_output.collect()
    path bams from perprocess_overview_bams_ch.collect{it[1]}
    path depth from samtools_depth_out.collect()
    path postdedup_stats from samtools_stats_postdedup_out.collect()
    path duplication_metric from mark_duplicate_metrics_ch.collect()
    path prededup_stats from samtools_stats_prededup_out.collect()
    path est_lib_comp from est_lib_comp_metrics_ch.collect()
    path input_tsv
    
    output:
    path "preprocessing_overview.tsv" into preprocess_summary_out_ch

    script:
    """
    python3 ${bin_dir}/pre_process_stats_overview.py \$PWD $input_tsv
    """
}
// END OF EVALUATION METRICS