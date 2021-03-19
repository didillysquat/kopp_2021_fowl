**A workflow for sWGS genotype likelihood calling and relatedness estimation (with no known variant sets)**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.07.1-brightgreen.svg)](https://www.nextflow.io/)
[![BINAC compatible](https://img.shields.io/badge/BINAC-compatible-brightgreen.svg)](https://www.nextflow.io/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**WARNING: This workflow and documentation are still in development**

- [Introduction](#introduction)
- [Quick Start](#quick-start)
  - [preprocess.nf quick start](#preprocessnf-quick-start)
  - [bqsr.nf quick start](#bqsrnf-quick-start)
  - [genotype.nf quick start](#genotypenf-quick-start)
- [Workflow summary](#workflow-summary)
  - [preprocess.nf summary](#preprocessnf-summary)
    - [preprocess.nf default pipeline](#preprocessnf-default-pipeline)
    - [preprocess.nf alternative options](#preprocessnf-alternative-options)
  - [bqsr.nf summary](#bqsrnf-summary)
    - [bqsr.nf default pipeline (without VCF provided as input)](#bqsrnf-default-pipeline-without-vcf-provided-as-input)
    - [bqsr.nf default pipeline (with VCF provided as input)](#bqsrnf-default-pipeline-with-vcf-provided-as-input)
  - [genotype.nf summary](#genotypenf-summary)
    - [genotype.nf default pipeline (Currently only implemented with GATK)](#genotypenf-default-pipeline-currently-only-implemented-with-gatk)
- [Documentation](#documentation)
  - [General documentation](#general-documentation)
    - [Nextflow arguments](#nextflow-arguments)
    - [General pipeline arguments](#general-pipeline-arguments)
    - [General pipeline outputs](#general-pipeline-outputs)
  - [preprocess.nf documentation](#preprocessnf-documentation)
    - [preprocess.nf arguments](#preprocessnf-arguments)
    - [preprocess.nf outputs](#preprocessnf-outputs)
  - [bqsr.nf documentation](#bqsrnf-documentation)
    - [bqsr.nf arguments](#bqsrnf-arguments)
    - [bqsr.nf outputs](#bqsrnf-outputs)
  - [genotype.nf documentation](#genotypenf-documentation)
    - [genotype.nf arguments](#genotypenf-arguments)
    - [genotype.nf outputs](#genotypenf-outputs)

## Introduction

This workflow is being developed as part of a research effort out of the [Kopp Group](https://www.ab.mpg.de/personen/98288) in collaboration with the [Sequencing Analysis Core Facility (SequAna)](https://www.biologie.uni-konstanz.de/sequana/) at the [University of Konstanz](https://www.uni-konstanz.de/en/).

The research effort is focused on relatedness in wild Guineafowl but this workflow may be applied to any diploid Eukaryotic organism for which an appropriate reference genome is available.

In particular the workflow is designed for organims for which there are no known/high confidence variant resources (e.g. dbSNP).

This workflow is made up of three nextflow pipelines that perform:

1. Preprocessing of input fastq to bam files - preprocess.nf
2. Base Quality Score Recalibration - bqsr.nf
3. Genotyping - genotype.nf

As an example use of the workflow, to generate a set of variants using two sets of BQSR you would run preprocess.nf, genotype.nf, bqsr.nf, genotype.nf, bqsr.nf, genotype.nf. As part of running this set of pipelines, you would be able to assess the effect of each round of BQSR. If further rounds were deemed necessary, another set of bqsr.nf and genotype.nf could be run.

If no BQSR is required, the bqsr.nf pipeline can be ommitted from the workflow.

## Quick Start

Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html).

Install one of [Singularity](https://singularity.lbl.gov/install-linux) of [Docker](https://docs.docker.com/get-docker/).

Provide your chosen container management system at the profile flag e.g.:
`-profile singularity` or `-profile docker`.

### preprocess.nf quick start

```bash
nextflow run preprocess.nf -profile <docker|singularity> --ref </path/to/your/ref/assembly.fna|fasta|fa> --input_dir </path/to/dir/with/paired/fastq.gz/reads/> --input_tsv </path/to/tsv/with/readgroup/info/>
```

[example readgroup info tsv]()

### bqsr.nf quick start

```bash
nextflow run bqsr.nf -profile <docker|singularity> --ref </path/to/your/ref/assembly.fna|fasta|fa> --input_dir </path/to/dir/with/bam/and/bam.bai/> [--input_vcf </path/to/multi-sample/vcf>]
```

### genotype.nf quick start

```bash
nextflow run genotyping.nf -profile <docker|singularity> --ref </path/to/your/ref/assembly.fna|fasta|fa> --input_dir </path/to/dir/with/bam/and/bam.bai/> 
```

## Workflow summary

### preprocess.nf summary

#### **preprocess.nf default pipeline**

- Create reference genome indices for mapping [[bwa mem](https://github.com/lh3/bwa)]
- Read quality assessment before processing [[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)]
- Adapter removal and sequence quality control [[trimmomatic](https://github.com/timflutre/trimmomatic)]
- Read quality assessment after processing [[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)]
- Read mapping to reference assembly [[bwa mem](https://github.com/lh3/bwa), [samtools](https://github.com/samtools)]
- Read merging [[gatk MergeSamFiles](https://gatk.broadinstitute.org/hc/en-us/articles/360046788832-MergeSamFiles-Picard-)]
- Annotation with read group headers [[gatk AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/360041850851-AddOrReplaceReadGroups-Picard-)]
- Remove duplicate reads [[gatk MarkDuplicatesSpark](https://gatk.broadinstitute.org/hc/en-us/articles/360057438771-MarkDuplicatesSpark)]
- Collect read metrics [GC metrics: [gatk CollectGcBiasMetrics](), PCR bottleneck coefficient, sequencing coverage: [mosdepth](https://github.com/brentp/mosdepth)]

#### **preprocess.nf alternative options**

**Subsampling:** Subsampling with [seqtk](https://github.com/lh3/seqtk) may be enabled for testing/dev purposes with a specified read depth using:

`--subsample --subsample_depth <number_of_seqs[int]>`

By default subsample_depth is set to 10000. 

**Mapping**: May be performed using [NextGenMap](https://github.com/Cibiv/NextGenMap).
Enabled with `--mapping ngm`
**Note:** due to an [issue](https://github.com/Cibiv/NextGenMap/issues/54) in the [ngm docker image](https://hub.docker.com/r/philres/nextgenmap/tags?page=1&ordering=last_updated), ngm is currently run using conda. As such, conda must be installed.

### bqsr.nf summary

#### **bqsr.nf default pipeline (without VCF provided as input)**

- Create reference genome dictionary [[gatk CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/360046222771-CreateSequenceDictionary-Picard-)]
- Create per sample GVCF [[gatk HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360050814612-HaplotypeCaller)]
- Create per sample GenomicsDB on a per reference scaffold basis [[gatk GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport)]
- Perform joint genotyping on a per reference scaffold basis [[gatk GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360046224151-GenotypeGVCFs)]
- Create high-confidence variant set on a per reference scaffold basis [[gatk VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360045800332-VariantFiltration), [gatk SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360047216851-SelectVariants)]
- Gather per scaffold high-confidence vcfs into single vcf [[gatk GatherVcfs](https://gatk.broadinstitute.org/hc/en-us/articles/360046787092-GatherVcfs-Picard-)]
- Make BQSR tables on pre-calibration BAMs [[gatk BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360050815072-BaseRecalibrator)]
- Apply BQSR tables [[gatk ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360046222011-ApplyBQSR)]
- Make BQSR tables on post-calibration BAMs [[gatk BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360050815072-BaseRecalibrator)]
- Calculate calibration metrics [[gatk AnalyzeCovariates](https://gatk.broadinstitute.org/hc/en-us/articles/360047215811-AnalyzeCovariates)]

#### **bqsr.nf default pipeline (with VCF provided as input)**

Same as above, but with the first three steps replaced with:
- Split multi-sample VCF to create per scaffold VCF set [[gatk SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360047216851-SelectVariants)]

### genotype.nf summary

#### **genotype.nf default pipeline (Currently only implemented with GATK)**

- Create reference genome dictionary [[gatk CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/360046222771-CreateSequenceDictionary-Picard-)]
- Create per sample GVCF [[gatk HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360050814612-HaplotypeCaller)]
- Create per sample GenomicsDB on a per reference scaffold basis [[gatk GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport)]
- Perform joint genotyping on a per reference scaffold basis [[gatk GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360046224151-GenotypeGVCFs)]
- Gather per scaffold vcfs into single vcf [[gatk GatherVcfs](https://gatk.broadinstitute.org/hc/en-us/articles/360046787092-GatherVcfs-Picard-)]
- Collect genotyping metrics [[bcftools stats](https://github.com/samtools/samtools), [rtg vcfstats](https://github.com/RealTimeGenomics/rtg-tools)]

## Documentation

### General documentation

#### **Nextflow arguments**

Nexflow specific options may be passed at execution. These start with a single hyphen (rather than a double hyphen).
Please refer to the [nextflow CLI documentation](https://www.nextflow.io/docs/latest/cli.html) for an exhaustive list.

Either `-profile singularity` or `-profile docker` must be provided at run time.
Singularity is more commonly used on HPC systems, while Docker is popular with private and single-user servers.

Passing `-resume` will make use of the automatically generated cache files.
If not passed the pipeline will rerun all processes.

In general, arguments may be passed on the command line (with a single hyphen for nextflow arguments or a double hyphen for pipeline-specific arguments) or by using a [configuration file](https://www.nextflow.io/docs/latest/config.html). `-C` specifies a specific configuration file to use, overriding any defaults. By default, the nextflow.config configuration file will be used.

#### **General pipeline arguments**

`--output_dir`: For each pipeline an output directory may be specified. If not specificed, a default output directory for each pipeline is used and will be created in the current working directory.

By default, the pipeline will not execute if an output directory already exists.
You can supply the `--overwrite` flag at execution to bypass this behaviour.
Sub directories for each output type will be created automatically.

#### **General pipeline outputs**
In addition to all specific output files for the given pipelines,
intermediate files generated by the pipeline are stored in the
`work` directory generated by Nextflow. This is found in the directory from which the pipeline was executed.

Nextflow will also produce a .nextflow.log file that contains the
logs for the current run.

### preprocess.nf documentation
The preprocess pipeline takes paired end fastq.gz files as input and outputs a set of analysis-read .bam and .bam.bai files.

#### **preprocess.nf arguments**
[Required Aguments]

`--input_dir`: The full path to the directory containing the input fastq.gz paired files

`--ref`: The full path to the reference genome assembly. It should be decompressed and in fasta format.

`--input_tsv`: The full path to the tab delimited file containing read group information. [example readgroup info tsv]()

[Optional Arguments]
TODO add memory and CPU optimisation arguments.

`--subsample`: Subsample each sample to the number of reads specified with `--subsample_depth <depth>`.
This option is useful during testing and development.

`--remake_indices`: Remake reference genome indices even if the genome indices already exist.

#### **preprocess.nf outputs**

**fastqc_pre_trim**: fastqc results before adapter trimming and sequence quality control

**fastqc_post_trim**: fastqc results after adapter trimming and sequence quality control

**markduplicates_metrics**: Metrics file produced by gatk MarkDuplicatesSpark

**output_bams**: The analysis-ready BAMs output from the pipeline.

**collect_gc_bias_metrics**: GC bias metrics

**pcr_bottleneck_coefficient**: PCR bottleneck calculation results

**mosdepth_sequencing_coverage**: Sequence coverage / depth results

### bqsr.nf documentation
Takes a set of BAM files as input and outputs a set of BAM files after applying a round of BQSR to them.

As part of the the GATK genotype calling we perform several rounds of BQSR. This workflow is designed to work without a curated set of known variants (as would normally be supplied to BQSR). Rather, a high confidence list of variants is generated through hard filtering an initially generated list of variants. BQSR is applied to the original BAMs using this list of variants to generate a new set of BAM files. These are then either used as the final set of BAM files from which to call variants for downstream analysis, or as input to a further round of BQSR. To help determine the optimal number of BQSR rounds we produce a metric summarising the effect of BQSR per round at the end of each round. Additionally, summary statistics are produced at the end of the genotype.nf pipeline such as Ti/Tv ratios that can be used to judge the quality of/improvement in the called variants and whether subsequent rounds of BQSR could be beneficial.

The fact that the BQSR process is cyclical in nature (ie. generating variants, to generate BAMs, to generate variants, to generate BAMs...), is why the workflow has been split up into three sepearte pipelines rather than a single continuous pipeline. Feedback loops have minimal support in Nextflow.

#### **bqsr.nf arguments**

[Required Aguments]

`--input_dir`: The full path to the directory containing the set of .bam and .bam.bai files to perform BQSR on. If running this after preprocess.nf, this will be the `output_bams` directory. If running this after a previous round of bqsr, this will be the `gatk_bqsr_output_bams` directory.

`--ref`: The full path to the reference genome assembly. It should be decompressed and in fasta format.

[Optional Arguments]
TODO add memory and CPU optimisation arguments.
`--input_vcf`: The full path to a multi-sample vcf file containing the called variants for the samples in question. If running bqsr.nf after genotype.nf this will be the .vcf file in the `gatk_output_variants` directory. N.B. the corresponding .vcf.idx file must be in the same directory. If supplied, the first three steps of the BQSR pipeline are skipped saving a considerable ammount of computational time. By making use of of an input vcf file, we save computing variants twice when running genotype.nf followed by bqsr.nf.

#### **bqsr.nf outputs**

TODO include the generated vcfs as output as they represent such a computational investment.
**gatk_bqsr_output_bams**: the post BQSR output .bam and .bam.bai files.

**gatk_bqsr_analyze_covariates_metrics**: Metrics summarising the effects of running BQSR. Generated by running gatk AnalyzeCovariates. See [this GitHub issue](https://github.com/broadinstitute/gatk/issues/322) for further details.

### genotype.nf documentation
genotype.nf takes a set of .bam and .bai files as input and outputs a multi-sample .vcf variant file containing genotype likelihood scores.

#### **genotype.nf arguments**

[Required Aguments]

`--input_dir`: The full path to the directory containing the set of .bam and .bam.bai files to call variants from. If running this after preprocess.nf, this will be the `output_bams` directory. If running this after bqsr.nf this will be the `gatk_bqsr_output_bams` directory.

`--ref`: The full path to the reference genome assembly. It should be decompressed and in fasta format.

[Optional Arguments]
TODO add memory and CPU optimisation arguments.

#### **genotype.nf outputs**

**gatk_output_vcf_publishDir**: The multi-sample .vcf file containing the genotype likelihoods.

**gatk_vcf_stats_publishDir**: Summary statistics generated for the called genotypes. These files are especially useful in deciding the number of BQSR rounds to incorporate into the analysis.
