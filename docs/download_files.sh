#!/usr/bin/env bash

# Download example files to test and validate hpexome

# Download BAM and index files from Illumina Platinum Genomes
wget \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12877_S1.bam \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12877_S1.bam.bai \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12878_S1.bam \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12878_S1.bam.bai \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12889_S1.bam \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12889_S1.bam.bai \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12890_S1.bam \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12890_S1.bam.bai \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12891_S1.bam \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12891_S1.bam.bai \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12892_S1.bam \
    https://storage.googleapis.com/genomics-public-data/platinum-genomes/bam/NA12892_S1.bam.bai

# Download known indels and sites VCF and index files
wget \
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

gunzip *

# Download reference genome files
wget \
    https://storage.googleapis.com/gatk-legacy-bundles/hg19/ucsc.hg19.dict \
    https://storage.googleapis.com/gatk-legacy-bundles/hg19/ucsc.hg19.fasta \
    https://storage.googleapis.com/gatk-legacy-bundles/hg19/ucsc.hg19.fasta.fai