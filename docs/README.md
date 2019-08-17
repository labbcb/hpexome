## hpexome: An automated workflow for processing whole-exome sequencing data

Whole-exome sequencing (WES) consists in the capture and sequence of all exons of protein coding genes in the human genome. It has been used in clinical applications like disease diagnosis. [Genome Analysis Toolkit (GATK)](https://software.broadinstitute.org/gatk/) is a programming framework for development of NGS data analysis tools. GATK provides tools for depth of coverage calculation analysis and single nucleotide polymorphism calling. These tools were used by large-scale sequencing projects such as 1000 Genomes Project and The Cancer Genome Atlas.

[GATK-Queue](http://gatkforums.broadinstitute.org/gatk/discussion/1306/overview-of-queue) is a command-line scripting framework for defining multi-stage genomic analysis pipelines combined with an execution manager that runs those pipelines from end-to-end. It has been used to implement workflows for WES data analysis. Our workflow for processing WES data was based on good practices described by authors of GATK software (Van der Auwera et al. 2002). It takes one or more BAM files, reference genome file, known SNP and indels databases, and outputs one or more VCF files.

Our workflow for processing WES data was based on good practices described by authors of GATK software. It takes one or more BAM files, reference genome file, known SNP and indels databases, and outputs one or more VCF files.

![Workflow for processing WES data](hpexome_workflow.png)

We used GATK-Queue workflow management to implement the workflow. It consists in a QScript file written in Scala programming language. The input parameters from every walker are set to default values but user can change theses values through the R interface.

## Installation

We provide a Python wrapper as command line tool for executing hpexome script. If not found this tool will download the latest version of Queue. To install our wrapper you must have installed Python 3 and pip.

``` bash
pip install --user git+https://github.com/labbcb/hpexome.git
```

## Usage

``` bash
hpexome [OPTIONS] [DESTINATION]
```

`DESTINATION` Sets the directory in which the outputs will be saved. If not set, the outputs will be saved in the directory in which the process is running.

OPTIONS:

__Required arguments__

- `-I, --bam` One or more sequence alignment files in BAM format _or_ directories containing `*.bam` files.
- `-R, --genome` Inputs a reference genetic sequence.
- `--dbsnp` Inputs the dbSNIP file, in which the IDs that are used to populate the ID column of the output.
- `--sites` Inputs a database of known polymorphic sites to skip over in the recalibration algorithm to be used in the BaseRecalibrator.

__Optional arguments__

- `--indels` Inputs the VCF file with known indels to be used with the Realigner Target Creator and the Indel Realigner.
- `-L, --intervals`Directs the GATK engine in restricting processing to genomic intervals. It is extremely useful when running the `scattercount` option with the Haplotype Caller, since it can prevent inaccuracies in its outputs.
- `--unified_vcf` If included, this parameter will result in the results of all inputs in the Haplotype Caller be joined into a single output.If not included, each input will generate its own output. Disabled by default.
- `-O, --output_file_name` Output file name for unified VCF. Default is `unified.vcf`.
- `--min_prunning` Sets the minimum support to not prune paths in the graph in the Haplotype Caller. If this is not specified, the program will use the default value, which is `2`.
- `-stand_call_conf` Sets the minimum phred-scaled confidence threshold at which variants should be called in the Haplotype Caller. If this not specified, the program will use the default value, which is __30__.

__Performance-specific arguments__

- ` --job_runner` Job executor engine (eg. Lsf706, Grid, PbsEngine).
- `-nt, --num_data_threads` Controls the number of data consecutive threads sent to the processor that are used in the parallelization process. It is used in the Realigner Target Creator, and may not be used together with the scattercount option. If not set, the walker will run in serial.
- `-nct, --num_threads_per_data_thread` Controls the number of CPU threads allocated to each data thread. It is used with the Base Recalibrator and the Print Reads, and may not be used together with the `scattercount` option. If not set, the walkers will run in serial.
- `--scatter_count` Controls the number of parts in which the genetic sequence will be divided when sent to be parallelized by the Job executor engine. It  is used in all walkers. It must be used with the `-jobRuner`  option, or else it will not use the GridEngine and the process will be run in serial.

__System path to required software__

- `--java_path` Path to java. Use this to pass JVM-specific arguments. Default is `java`.
- `--queue_path` Path to Queue jar file. Default is `Queue.jar`.

# Examples

Using local computing environment

``` bash
hpexome \
    --bam NA12877_S1.bam \
    --bam NA12878_S1.bam \
    --genome ucsc.hg19.fasta  \
    --dbsnp dbsnp_138.hg19.vcf \
    --indels Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    --indels 1000G_phase1.indels.hg19.sites.vcf \
    --sites 1000G_phase1.snps.high_confidence.hg19.sites.vcf \
    --sites 1000G_omni2.5.hg19.sites.vcf \
    --unified_vcf \
    --num_data_threads 4 \
    --num_threads_per_data_thread 4 \
    results
```

Using cluster environment (PBS/Torque)

``` bash
hpexome \
    --bam NA12877_S1.bam \
    --bam NA12878_S1.bam \
    --genome ucsc.hg19.fasta  \
    --dbsnp dbsnp_138.hg19.vcf \
    --indels Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    --indels 1000G_phase1.indels.hg19.sites.vcf \
    --sites 1000G_phase1.snps.high_confidence.hg19.sites.vcf \
    --sites 1000G_omni2.5.hg19.sites.vcf \
    --unified_vcf \
    --scatter_count 4 \
    --job_runner GridEngine \
    results
```

## References

__From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline__ Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M, 2013 _CURRENT PROTOCOLS IN BIOINFORMATICS 43:11.10.1-11.10.33_. [Article](http://dx.doi.org/10.1002/0471250953.bi1110s43)
