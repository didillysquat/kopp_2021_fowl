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
This workflow takes two vcf files as input (GATK- and bcftools-based),
computes the variants in common between the two vcfs using bcftools isec,
excludes the mitochondrial scaffold (vcftools),
and thins out the variants (vcftools) using filters of --maf 0.025, --max-missing 0.75, --remove-filtered-all --remove-indels.
Relatedness is then calculated using three different pieces of software:
lcmlkin: https://github.com/didillysquat/maximum-likelihood-relatedness-estimation
read: https://bitbucket.org/tguenther/read/src/master/
NgsRelate: https://github.com/ANGSD/NgsRelate
*/

gatk_vcfgz = file(params.gatk_vcfgz)
gatk_vcfgz_tbi = file("${gatk_vcfgz}.tbi")
bcftools_vcfgz = file(params.bcftools_vcfgz)
bcftools_vcfgz_tbi = file("${bcftools_vcfgz}.tbi")
bin_dir = "${workflow.launchDir}/bin"
params.output_dir = "${workflow.launchDir}/outputs/relatedness"
pca_publish_dir = [params.output_dir, "PCA"].join(File.separator)
read_publish_dir = [params.output_dir, "read"].join(File.separator)
lcmlkin_publish_dir = [params.output_dir, "lcmlkin"].join(File.separator)
ngsrelate_publish_dir = [params.output_dir, "ngsrelate"].join(File.separator)
thinned_vcf_publish_dir = [params.output_dir, "thinned_vcf"].join(File.separator)

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

if (!(params.isec_threads)){
    params.isec_threads = 1
}

if (!(params.lcmlkin_threads)){
    params.lcmlkin_threads = 1
}

if (!(params.ngsrelate_threads)){
    params.ngsrelate_threads = 1
}

if (params.bam_dir_for_PCA){
    def path_to_check = new File(params.bam_dir_for_PCA);
    if (!path_to_check.exists()) {
        throw new Exception("""The bam_dir_for_PCA directory ${path_to_check} doesn't exists.\n
        Please provide the full path to the bam directory to use \n""")
    }
    do_pca = true
    
    // We will generate the group name by removing all parts of the file name that are found in common between all samples.
    pca_bams_in_ch = Channel.fromFilePairs("${params.bam_dir_for_PCA}/*.bam{,.bai}").toList().flatMap{
        list_of_bams -> 
        // println("This is the list_of_bams: ${list_of_bams}")
        def group_list = [];
        // For each of the bam sets, pull out the current group name and put it into the group_list
        list_of_bams.eachWithIndex{ bam_set, i -> group_list << bam_set[0];}
        // println("This is the group_list: ${group_list}" )
        // Then work backwards on a per character basis checking to see if the character is found in common
        // for all samples. If it is then this will be a character to discard. Keep going until
        // we find a character that is not in common across all names and make a note of the index position.
        // We will then return the group ids with this number of characters discarded
        def cut_index = 0;
        def exit = false;
        for (i = 1; !exit; i++) {
            // For each character starting from the last
            start_char = group_list[0].charAt(group_list[0].length() - i)
            // println("Start char is currently: ${start_char}")
            group_list.eachWithIndex{
                group_name, j -> 
                    if (group_name.charAt(group_name.length() - i) != start_char){
                        // Then j is the index to start cutting from
                        cut_index = -1 * i;
                        // println("Look we got here and the cut index is: ${cut_index}")
                        exit = true;
                    }
                }
        }
        // Now that we have the cut index we can cut each of the group names
        // At this point we want to create a new list and return it
        def output_list = [];
        list_of_bams.eachWithIndex{bam_set, i -> output_list << [bam_set[0][0..cut_index], bam_set[1] ];}
        return output_list;
    }
}else{
    print("No bam_dir_for_PCA provided. PCA will not be conducted.")
    do_pca = false
}

process isec{
    container "halllab/bcftools:v1.9"
    cpus params.isec_threads

    input:
    file gatk_vcfgz
    file gatk_vcfgz_tbi
    file bcftools_vcfgz
    file bcftools_vcfgz_tbi

    output:
    path("relatedness.isec.0002.vcf") into exclude_mito_scaff_ch

    script:
    """
    mkdir isec
    bcftools isec $bcftools_vcfgz $gatk_vcfgz -p isec --threads ${task.cpus}
    mv isec/0002.vcf relatedness.isec.0002.vcf
    """
    
}

// If a mitochondrial scaffold name is given, exclude this before thinning
if (params.mito_scaff){
    process exclude_mito_scaff{
        container  "biocontainers/vcftools:v0.1.16-1-deb_cv1"
        if (workflow.containerEngine == 'docker'){
            containerOptions '-u $(id -u):$(id -g)'
        }

        input:
        path(isec_vcf) from exclude_mito_scaff_ch

        output:
        path "relatedness.isec.0002.exMito.vcf" into thin_ch

        script:
        """
        vcftools --not-chr ${params.mito_scaff} --recode --recode-INFO-all --stdout --gzvcf $isec_vcf > relatedness.isec.0002.exMito.vcf
        """
    }
}else{
    exclude_mito_scaff.into{thin_ch}
}

process thin{
    container "biocontainers/vcftools:v0.1.16-1-deb_cv1"
    publishDir thinned_vcf_publish_dir, mode: "copy"
    if (workflow.containerEngine == 'docker'){
            containerOptions '-u $(id -u):$(id -g)'
        }

    input:
    path vcf_to_thin from thin_ch

    output:
    path "relatedness.isec.0002.exMito.thinned.vcf" into lcmlkin_ch,read_ch,ngsrelate_ch

    script:
    """
    vcftools --vcf $vcf_to_thin --remove-filtered-all --remove-indels --maf 0.025 --recode --recode-INFO-all --stdout --max-missing 0.75 > relatedness.isec.0002.exMito.thinned.vcf
    """
}

if (do_pca){
    
    // Generate the beagle formated genotype likelihood files for input to pcangsd
    process generate_PCAngsd_input{
        tag "generate_PCAngsd_input"
        container 'autamus/angsd:latest'
        cpus 4

        input:
        file bams from pca_bams_in_ch.collect(){it[1]}

        output:
        file "genolike.beagle.gz" into pcangsd_in_ch
        file "sample.names.order.txt" into pc_sample_order_ch

        shell:
        '''
        find ~+ -name '*.bam' > bam_list.full.txt
        awk -F'/' '{print $NF}' bam_list.full.txt > bam_list.short.txt
        awk -f !{bin_dir}/unique_names.awk bam_list.short.txt > sample.names.order.txt
        angsd -GL 2 -out genolike -P I{task.cpus} -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam bam_list.full.txt
        '''
    }

    // Generate the covariance matrix
    process generate_PCAngsd_cov_matrix{
        tag "generate_PCAngsd_cov_matrix"
        container 'didillysquat/pcangsd:latest'
        cpus 4

        input:
        file genolike from pcangsd_in_ch

        output:
        file "pcangsd.out.cov" into generate_pcs_ch

        script:
        """
        pcangsd.py -beagle genolike.beagle.gz -out pcangsd.out -threads ${task.cpus}
        """
    }

    // Eigen decompose the covariance matrix to generate the PCs
    // Also annotate with the sample names so user knows the order
    process generate_PCs{
        tag "generate_PCs"
        container 'didillysquat/r_multipurpose:latest'
        publishDir pca_publish_dir

        input:
        file cov_file from generate_pcs_ch
        file pc_sample_order from pc_sample_order_ch

        output:
        tuple file("PCs.txt"), file("PC.values.txt"), file("PC.plot.png") into generate_pcs_out_ch

        script:
        """
        Rscript ${bin_dir}/compute_pcs.r \$PWD $pc_sample_order $cov_file
        """
    }
}


process read_make_tped_tfam{
    container "biocontainers/vcftools:v0.1.16-1-deb_cv1"
    if (workflow.containerEngine == 'docker'){
            containerOptions '-u $(id -u):$(id -g)'
        }

    input:
    path read_vcf from read_ch

    output:
    tuple path("relatedness.isec.0002.exMito.thinned.tped"), path("relatedness.isec.0002.exMito.thinned.tfam") into read_run_ch

    script:
    """
    vcftools --vcf $read_vcf --plink-tped --out relatedness.isec.0002.exMito.thinned
    """
}

process read_run{
    container "didillysquat/read:latest"
    publishDir read_publish_dir, mode: "copy"

    input:
    tuple path(tped), path(tfam) from read_run_ch

    output:
    tuple path("READ_results"), path("READ_results_plot.pdf"), path("READ_output_ordered"), path("Read_intermediate_output"), path("meansP0_AncientDNA_normalized") into read_out_ch

    script:
    base_name = tped.getName().replaceAll(".tped", "")
    """
    READ.py $base_name
    """
}

process lcmlkin{
    container "didillysquat/lcmlkin:latest"
    publishDir lcmlkin_publish_dir, mode: "copy"
    cpus params.lcmlkin_threads

    input:
    path vcf_in from lcmlkin_ch

    output:
    tuple path("relatedness.isec.0002.exMito.thinned.lcmlkin.results"), path("relatedness.isec.0002.exMito.thinned.lcmlkin.results.log") into lcmlkin_out_ch

    script:
    """
    lcmlkin -g all --likelihoodFormat phred -o relatedness.isec.0002.exMito.thinned.lcmlkin.results -i $vcf_in -t ${task.cpus}
    """
}

process ngsrelate_threads{
    container "didillysquat/ngsrelate:latest"
    publishDir ngsrelate_publish_dir, mode: "copy"
    cpus params.ngsrelate_threads

    input:
    path vcf_in from ngsrelate_ch

    output:
    path("relatedness.isec.0002.exMito.thinned.ngsrelate.results") into ngsrelate_threads_out_ch

    script:
    """
    awk '/#CHROM/ {for(i=10; i<=NF; i++){print \$i}}' relatedness.isec.0002.exMito.thinned.vcf > sample_names.txt
    ngsRelate -h $vcf_in -O relatedness.isec.0002.exMito.thinned.ngsrelate.results -c 1 -p ${task.cpus} -z sample_names.txt
    """
}