#!/usr/bin/env nextflow

params.raw_reads_dir = "/home/humebc/projects/20210125_kopp_guinea_fowl/raw_seq_files"
params.ref_assembly_path = "/home/humebc/projects/20210125_kopp_guinea_fowl/nf_pipeline_kopp/nummel_ref_assembly/GCF_002078875.1_NumMel1.0_genomic.fna.gz"

params.path_to_read_group_lists = "/home/humebc/projects/20210125_kopp_guinea_fowl/raw_seq_files/file-list-run-1,/home/humebc/projects/20210125_kopp_guinea_fowl/raw_seq_files/file-list-run-2"


bin_dir = "${workflow.launchDir}/bin"
launch_dir = "${workflow.launchDir}"
tru_seq_pe_fasta_path = "${workflow.launchDir}/TruSeq3-PE.fa"

params.subsample = true
params.subsample_depth = 10000

// TODO make a dictionary that maps the pair_id to the readgroup info we will need to add readgroups using
// picard AddOrReplaceReadGroups
// This will return a map with key that is the pair_id (taking into account whether we are subsampling
// and therefore have a modified paid_id) and the value is the string that will be passed to the actual
// picard AddOrReplaceReadGroups command e.g. RGID=xxx RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=xxx
// NB it is critical to convert the key for a map to a string using .toString() if using e.g. "${my_var}"
read_group_map = { 
    String[] file_lists = params.path_to_read_group_lists.split(",");
    def emptyMap = [:]
    file_lists.eachWithIndex{ file_list, i -> 
            // i will be the RGID value
            // j will be the RGSM value
        def j = 0
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
                            emptyMap.put("${pair_id}".toString(), "RGID=${i+1} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${j}")
                        }                            
                    }else{
                        throw new Exception("An error has occured when making the readgroup dictionary\nNo match was found for ${file_name}")
                    }
                }
            }
        }
    }
    return emptyMap     
}()


// def bob = "mpg_L16952-1_W1446_S28_sub_10000"
// println(read_group_map["${bob}"])

// println("This is 'mpg_L16952-1_W1446_S28_sub_10000': ${read_group_map[mpg_L16952-1_W1446_S28_sub_10000]}")

// Publish dirs
if (params.subsample){
    params.output_dir = "${workflow.launchDir}/outputs_sub_sampled/"
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim_sub_${params.subsample_depth}"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim_sub_${params.subsample_depth}"].join(File.separator)
}else{
    params.output_dir = "${workflow.launchDir}/outputs"
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim"].join(File.separator)
}

// CPUs
params.trimmomatic_threads = 50
params.nextgenmap_threads = 40
params.markduplicates_cpus = 20


// /* 
// If subsample, create a channel that will pass into the subsampling process and then pass the subsampled
// files into the trimming and fastqc
// Modify the pair name to indicate subsampled files. This way it will carry through the remainder of the anlysis.
// if not subsample, then create channel directly from the raw sequencing files
// */
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


/* 
    Trim reads with trimmomatic
*/
/*  Call used by Till:
        trimmomatic PE -threads 8 \ 
		${infile} ${base}R2_001.fastq.gz \
                ${base}R1_001_paired.fastq.gz ${base}R1_001_unpaired.fastq.gz \
                ${base}R2_001_paired.fastq.gz ${base}R2_001_unpaired.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:70 LEADING:10 HEADCROP:5 TRAILING:10 \
		-summary ${base}trimSum.txt
*/
// TODO check to see if we need to have a higher HEADCROP than Till used.
// We need to pick up the paired and unpaired reads
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
    // tuple val(pair_id), file("${pair_id}*{1,2}{U,P}fq.gz") into ch_fastqc_post_trim

	script:
	outbase = fastqs[0].getName().replaceAll('.fastq.gz', '.trimmed.fq.gz')
	"""
	trimmomatic PE -threads ${params.trimmomatic_threads} -basein ${fastqs[0]} \\
		-baseout $outbase \\
		ILLUMINACLIP:${tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
		SLIDINGWINDOW:4:20 MINLEN:70 LEADING:10 HEADCROP:5 TRAILING:10 
	"""
}

// ch_fastqc_post_trim_paired.join(ch_fastqc_post_trim_unpaired).view()

// We will look specifically for the P1 P2 U1 and U2 files
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

// We don't use the genome variable that comes from the index_ch
// This input is only there to ensure that the indexing has been performed before the mapping
// TODO possibly remove the all_reads. I'm not sure if we'll need this.
process nextgenmap_mapping{
    tag "${pair_id}"
    conda 'envs/ngm.yaml'
    cpus params.nextgenmap_threads

    input:
    path genome from index_ch
	tuple val(pair_id), file(paired), file(unpaired) from ch_ngm_paired.join(ch_ngm_unpaired)

    output:
    tuple val(pair_id), file("${pair_id}_P_mapped.sam"), file("${pair_id}_U1_mapped.sam"), file("${pair_id}_U2_mapped.sam"), file("${pair_id}_allreads_mapped.sam") into add_read_group_headers_ch

    script:
    """
    ngm -t ${task.cpus} -r ${params.ref_assembly_path} -1 ${paired[0]} -2 ${paired[1]} -o ${pair_id}_P_mapped.sam;
    ngm -t ${task.cpus} -r ${params.ref_assembly_path} -q ${unpaired[0]} -o ${pair_id}_U1_mapped.sam;
    ngm -t ${task.cpus} -r ${params.ref_assembly_path} -q ${unpaired[1]} -o ${pair_id}_U2_mapped.sam;
    cat *mapped.sam > ${pair_id}_allreads_mapped.sam
    """
}

// TODO here we need to add the readgroup headers
process add_read_group_headers{
    tag "${pair_id}"
    container 'broadinstitute/gatk:latest'
    maxForks 1

    input:
    tuple val(pair_id), file(paired), file(unpaired_one), file(unpaired_two), file(all_reads) from add_read_group_headers_ch

    output:
    tuple val(pair_id), file("${pair_id}.paired.mapped.headers.sam"), file("${pair_id}.unpaired.1.mapped.headers.sam"), file("${pair_id}.unpaired.2.mapped.headers.sam") into mark_duplicates_ch

    script:
    """
    gatk AddOrReplaceReadGroups I=${paired} O=${pair_id}.paired.mapped.headers.sam ${read_group_map[(pair_id)]}
    gatk AddOrReplaceReadGroups I=${unpaired_one} O=${pair_id}.unpaired.1.mapped.headers.sam ${read_group_map[(pair_id)]}
    gatk AddOrReplaceReadGroups I=${unpaired_two} O=${pair_id}.unpaired.2.mapped.headers.sam ${read_group_map[(pair_id)]}
    """
}

// // TODO possibly remove the all_reads. I'm not sure if we'll need this.
// // containerOptions '-u $(id -u):$(id -g)'
// process markduplicates_spark{
//     tag "${pair_id}"
//     container 'broadinstitute/gatk:latest'
//     cpus params.markduplicates_cpus

//     input:
//     tuple val(pair_id), file(paired), file(unpaired_fwd), file(unpaired_rev), file(all_reads) from mark_duplicates_ch

//     output:
//     tuple val(pair_id), file("${pair_id}.paired.deduplicated.sorted.bam"), file("${pair_id}.unpaired_fwd.deduplicated.sorted.bam"), file("${pair_id}.unpaired_rev.deduplicated.sorted.bam")
//     tuple file("${pair_id}.paired.deduplicated.sorted.metrics.txt"), file("${pair_id}.unpaired_fwd.deduplicated.sorted.metrics.txt"), file("${pair_id}.unpaired_rev.deduplicated.sorted.metrics.txt")

//     script:
//     """
//     gatk MarkDuplicatesSpark --remove-sequencing-duplicates --conf 'spark.executor.cores=${task.cpus}' -I ${paired} -O ${pair_id}.paired.deduplicated.sorted.bam -M ${pair_id}.paired.deduplicated.sorted.metrics.txt
//     gatk MarkDuplicatesSpark --remove-sequencing-duplicates --conf 'spark.executor.cores=${task.cpus}' -I ${unpaired_fwd} -O ${pair_id}.unpaired_fwd.deduplicated.sorted.bam -M ${pair_id}.unpaired_fwd.deduplicated.sorted.metrics.txt
//     gatk MarkDuplicatesSpark --remove-sequencing-duplicates --conf 'spark.executor.cores=${task.cpus}' -I ${unpaired_rev} -O ${pair_id}.unpaired_rev.deduplicated.sorted.bam -M ${pair_id}.unpaired_rev.deduplicated.sorted.metrics.txt
//     """
// }