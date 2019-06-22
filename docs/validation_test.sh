#!/usr/bin/env bash

CSV_FILE=hpexome_timings.csv
ALN_DIR=~/bioinf/aln
REF_DIR=~/bioinf/ref/hg19
RESULT_DIR=~/bioinf/res

hpexome \
    --bam ${ALN_DIR}/NA12877_S1.bam \
    --genome ${REF_DIR}/ucsc.hg19.fasta \
    --dbsnp ${REF_DIR}/dbsnp_138.hg19.vcf \
    --indels ${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    --indels ${REF_DIR}/1000G_phase1.indels.hg19.sites.vcf \
    --sites ${REF_DIR}/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
    --sites ${REF_DIR}/1000G_omni2.5.hg19.sites.vcf \
    --scatter_count 16 \
    --job_runner PbsEngine \
    ${RESULT_DIR}/hpexome_validation