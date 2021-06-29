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
params.output_dir = "${workflow.launchDir}/outputs/relatedness"
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
    tuple path("READ_results"), path("READ_results_plot.pdf") into read_out_ch

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
    ngsRelate -h $vcf_in -O relatedness.isec.0002.exMito.thinned.ngsrelate.results -c 1 -p ${task.cpus}
    """
}