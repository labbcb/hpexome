# HPexome

**HPexome** only requires Python 3 and Java 8 to run, optionally DMRAA-supported batch processing system such as SGE.
However, to create input data it is required to align raw sequencing reads to reference genome and sort those reads by coordinate.
The required software are: BWA to align reads, amtools to convert output to BAM, and Picard to sort reads and fix tags.

```bash
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

# HPexome
python3 -m venv venv
source venv/bin/activate
pip install hpexome
```

## Performance test

Using [NA12878 exome data from 1000 Genomes phase 3 release](https://www.internationalgenome.org/data-portal/sample/NA12878).

Download and generate required input files.

```bash
# Download FASTQ files
wget --quiet -i fastq_files.txt

# Download and unzip reference files
wget --quiet -i reference_files.txt
gunzip dbsnp_138.hg19.vcf.gz
    dbsnp_138.hg19.vcf.idx.gz
    Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
    Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz
    1000G_phase1.indels.hg19.sites.vcf.gz
    1000G_phase1.indels.hg19.sites.vcf.idx.gz
    1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
    1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx.gz
    1000G_omni2.5.hg19.sites.vcf.gz
    1000G_omni2.5.hg19.sites.vcf.idx.gz

# Build BWA genome index and align reads
bwa-0.7.17/bwa index ucsc.hg19.fasta
bwa-0.7.17/bwa mem \
    -K 100000000 -t 16 -Y ucsc.hg19.fasta \
    SRR098401_1.filt.fastq.gz SRR098401_2.filt.fastq.gz \
    | samtools-1.10/samtools view -1 - > NA12878.bam

# Sort aligned reads by genomic coordinates
java -jar picard.jar SortSam \
    INPUT=NA12878.bam \
    OUTPUT=NA12878.sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

# Fix ReadGroup tags
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

File tree.

    hpexome-paper/
    ├── 1000G_omni2.5.hg19.sites.vcf
    ├── 1000G_omni2.5.hg19.sites.vcf.idx
    ├── 1000G_phase1.indels.hg19.sites.vcf
    ├── 1000G_phase1.indels.hg19.sites.vcf.idx
    ├── 1000G_phase1.snps.high_confidence.hg19.sites.vcf
    ├── 1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx
    ├── bwa-0.7.17
    ├── bwa-0.7.17.tar.bz2
    ├── dbsnp_138.hg19.vcf
    ├── dbsnp_138.hg19.vcf.idx
    ├── fastq_files.txt
    ├── Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
    ├── Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx
    ├── NA12878.bam
    ├── NA12878.sorted.bai
    ├── NA12878.sorted.bam
    ├── NA12878.sorted.rgfix.bai
    ├── NA12878.sorted.rgfix.bam
    ├── performance.Rmd
    ├── picard.jar
    ├── README.md
    ├── reference_files.txt
    ├── run_performance_tests.sh
    ├── samtools-1.10
    ├── samtools-1.10.tar.bz2
    ├── SRR098401_1.filt.fastq.gz
    ├── SRR098401_2.filt.fastq.gz
    ├── ucsc.hg19.dict
    ├── ucsc.hg19.fasta
    ├── ucsc.hg19.fasta.amb
    ├── ucsc.hg19.fasta.ann
    ├── ucsc.hg19.fasta.bwt
    ├── ucsc.hg19.fasta.fai
    ├── ucsc.hg19.fasta.pac
    ├── ucsc.hg19.fasta.sa
    ├── validation.Rmd
    └── venv

Run performance test script.

```bash
export SGE_ROOT=/var/lib/gridengine
bash run_performance_tests.sh
```

The script generated [hpexome_timings.csv](hpexome_timings.csv) file.

Report is available at [performance.md](performance.md).

## Validation test

```bash
export SGE_ROOT=/var/lib/gridengine
hpexome \
    --bam NA12877_S1.bam \
    --genome ucsc.hg19.fasta \
    --dbsnp dbsnp_138.hg19.vcf \
    --indels ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz \
    --indels ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf.gz \
    --sites ALL.wgs.dbsnp.build135.snps.sites.vcf.gz \
    --stand_call_conf 40 \
    --stand_emit_conf 10 \
    --min_pruning 3 \
    --scatter_count 16 \
    --job_runner GridEngine \
    hpexome_validation
```

Validation report is available at [validation.md](validation.md).
