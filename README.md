# An automated tool for processing whole-exome sequencing data

Whole-exome sequencing has been widely used in clinical applications for the identification of the genetic causes of several diseases.
_HPexome_ automates many data processing tasks for exome-sequencing data analysis of large-scale cohorts.
Given ready-analysis alignment files it is capable of breaking input data into small genomic regions to efficiently process in parallel on cluster-computing environments.
It relies on Queue workflow execution engine and GATK variant calling tool and its best practices to output high-confident unified variant calling file.
Our workflow is shipped as Python command line tool making it easy to install and use.

``` bash
hpexome [OPTIONS] [DESTINATION]
```

`OPTIONS`

__Required arguments__

- `-I, --bam` One or more sequence alignment files in BAM format _or_ directories containing `*.bam` files.
- `-R, --genome` Reference genome in single FASTA file.
- `--dbsnp` dbSNP file in VCF format. See [dbSNP FTP](ftp://ftp.ncbi.nih.gov/snp/).
- `--sites` VCF files containing known polymorphic sites to skip over in the recalibration algorithm.

__Optional arguments__

- `--indels` Inputs the VCF file with known insertions and deletions (indels) to be used.
- `-L, --intervals` One or more genomic intervals over which to operate.
- `--unified_vcf` Unify VCF files into a single one.
- `-O, --output_file_name` Output file name for unified VCF. Default is `unified.vcf`.
- `--min_prunning` Minimum support to not prune paths in the graph. Default value is `2`.
- `-stand_call_conf` Minimum phred-scaled confidence threshold at which variants should be called. Default is `30`.

__Performance-specific arguments__

- `-nt, --num_data_threads` Controls the number of data consecutive threads sent to the processor that are used in the parallelization process. It is used in the Realigner Target Creator, and may not be used together with the scattercount option. If not set, the walker will run in serial.
- `-nct, --num_threads_per_data_thread` Controls the number of CPU threads allocated to each data thread. It is used with the Base Recalibrator and the Print Reads, and may not be used together with the `scattercount` option. If not set, the walkers will run in serial.
- ` --job_runner` Job executor engine (eg. Lsf706, Grid, PbsEngine).
- `--scatter_count` Controls the number of parts in which the genetic sequence will be divided when sent to be parallelized by the Job executor engine. It  is used in all walkers. It must be used with the `-jobRuner`  option, or else it will not use the GridEngine and the process will be run in serial.

__System path to required software__

- `--java_path` Path to java. Use this to pass JVM-specific arguments. Default is `java`.

`DESTINATION` Sets the directory in which the outputs will be saved. If not set, the outputs will be saved in the directory in which the process is running.

## Reproducible example

In this example we will download and process a [whole-exome sequence sample](https://www.internationalgenome.org/data-portal/sample/HG00114) from the 1000 Genomes Project and required reference files, as well as required software.
Let's create a directory to write all files.

```bash
mkdir hpexome
cd hpexome
```

**HPexome** only requires Python 3 and Java 8 to run, optionally DMRAA-supported batch processing system such as SGE.
However, to create input data it is required to align raw sequencing reads to reference genome and sort those reads by coordinate.
The required software are: BWA to align reads, amtools to convert output to BAM, and Picard to sort reads and fix tags.

```bash
# HPexome
pip3 install hpexome

# BWA
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar xf bwa-0.7.17.tar.bz2 
make -C bwa-0.7.17 

# Samtools
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar xf samtools-1.10.tar.bz2
make -C samtools-1.10

# Picard
wget https://github.com/broadinstitute/picard/releases/download/2.21.7/picard.jar
```

Download raw sequencing data as FASTQ files.

```bash
wget \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/sequence_read/SRR098401_1.filt.fastq.gz \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/sequence_read/SRR098401_2.filt.fastq.gz
```

Download required reference files.

```bash
wget \
    https://storage.googleapis.com/gatk-legacy-bundles/hg19/ucsc.hg19.fasta \
    https://storage.googleapis.com/gatk-legacy-bundles/hg19/ucsc.hg19.fasta.fai \
    https://storage.googleapis.com/gatk-legacy-bundles/hg19/ucsc.hg19.dict \
    ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz \
    ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz \
    ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \
    ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz \
    ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz \
    ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz \
    ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
    ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx.gz \
    ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz \
    ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.idx.gz

gunzip dbsnp_138.hg19.vcf.gz \
    dbsnp_138.hg19.vcf.idx.gz \
    Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \
    Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz \
    1000G_phase1.indels.hg19.sites.vcf.gz \
    1000G_phase1.indels.hg19.sites.vcf.idx.gz \
    1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
    1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx.gz \
    1000G_omni2.5.hg19.sites.vcf.gz \
    1000G_omni2.5.hg19.sites.vcf.idx.gz
```

Index reference genome.

```bash
bwa-0.7.17/bwa index ucsc.hg19.fasta
```

Align raw sequencing reads to the human reference genome.

```bash
bwa-0.7.17/bwa mem \
    -K 100000000 -t 16 -Y ucsc.hg19.fasta \
    SRR098401_1.filt.fastq.gz SRR098401_2.filt.fastq.gz \
    | samtools-1.10/samtools view -1 - > NA12878.bam
```

Sort aligned reads by genomic coordinates.

```bash
java -jar picard.jar SortSam \
    INPUT=NA12878.bam \
    OUTPUT=NA12878.sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true
```

Fix RG tags.

```bash
java -jar picard.jar AddOrReplaceReadGroups \
    I=NA12878.sorted.bam \
    O=NA12878.sorted.rgfix.bam \
    RGID=NA12878 \
    RGSM=NA12878 \
    RGLB=1kgenomes \
    RGPL=Illumina \
    PU=Unit1 \
    CREATE_INDEX=true
```

In some computing setups it will require to set `SGE_ROOT` environment variable.

```bash
export SGE_ROOT=/var/lib/gridengine
```

Run **HPexome**.

```bash
hpexome \
    --bam NA12878.sorted.rgfix.bam \
    --genome ucsc.hg19.fasta \
    --dbsnp dbsnp_138.hg19.vcf \
    --indels Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    --indels 1000G_phase1.indels.hg19.sites.vcf \
    --sites 1000G_phase1.snps.high_confidence.hg19.sites.vcf \
    --sites 1000G_omni2.5.hg19.sites.vcf \
    --scatter_count 16 \
    --job_runner GridEngine \
    result_files
```

It is expected the following files.

    result_files/
    ├── HPexome.jobreport.txt
    ├── NA12878.sorted.rgfix.HC.raw.vcf
    ├── NA12878.sorted.rgfix.HC.raw.vcf.idx
    ├── NA12878.sorted.rgfix.HC.raw.vcf.out
    ├── NA12878.sorted.rgfix.intervals
    ├── NA12878.sorted.rgfix.intervals.out
    ├── NA12878.sorted.rgfix.realn.bai
    ├── NA12878.sorted.rgfix.realn.bam
    ├── NA12878.sorted.rgfix.realn.bam.out
    ├── NA12878.sorted.rgfix.recal.bai
    ├── NA12878.sorted.rgfix.recal.bam
    ├── NA12878.sorted.rgfix.recal.bam.out
    ├── NA12878.sorted.rgfix.recal.cvs
    └── NA12878.sorted.rgfix.recal.cvs.out

See [/docs](/docs) folder for information about performance and validation tests.

## Development

Clone git repository from GitHub.

```bash
git clone https://github.com/labbcb/hpexome.git
cd hpexome
```

Create virtual environment and install development version.

```bash
python3 -m venv venv
source venv/bin/activate
pip install --requirement requirements.txt
```

Publish new hpexome version to Pypi.

```bash
pip install setuptools wheel twine
python setup.py sdist bdist_wheel
twine upload -u $PYPI_USER -p $PYPI_PASS dist/*
```
