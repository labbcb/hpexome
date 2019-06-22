#!/usr/bin/env bash

# Run hpexome workflow changing scatter count, repeat 5 times each
# Scatter count values: 1, 2, 4, 8, 16

CSV_FILE=hpexome_timings.csv
ALN_DIR=~/bioinf/aln
REF_DIR=~/bioinf/ref/hg19
RESULT_DIR=~/bioinf/res

echo "scatter_count, count, start, end" > ${CSV_FILE}
for scatter_count in 1 2 4 8 16 32 64
do
    for count in 1 2 3 4 5
    do
        start="$(date +%s)"
        hpexome \
            --bam ${ALN_DIR}/NA12877_S1.bam \
            --bam ${ALN_DIR}/NA12878_S1.bam \
            --bam ${ALN_DIR}/NA12889_S1.bam \
            --bam ${ALN_DIR}/NA12890_S1.bam \
            --bam ${ALN_DIR}/NA12891_S1.bam \
            --bam ${ALN_DIR}/NA12892_S1.bam \
            --genome ${REF_DIR}/ucsc.hg19.fasta \
            --dbsnp ${REF_DIR}/dbsnp_138.hg19.vcf \
            --indels ${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
            --indels ${REF_DIR}/1000G_phase1.indels.hg19.sites.vcf \
            --sites ${REF_DIR}/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
            --sites ${REF_DIR}/1000G_omni2.5.hg19.sites.vcf \
            --scatter_count ${scatter_count} \
            --job_runner PbsEngine \
            ${RESULT_DIR}/hpexome_${scatter_count}_${count}

        end="$(date +%s)"
        echo "${scatter_count}, ${count}, ${start}, ${end}" > ${CSV_FILE}
    done
done
