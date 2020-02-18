#!/usr/bin/env bash

# Run hpexome workflow changing scatter count, repeat 5 times each
# Scatter count values: 1, 2, 4, 8, 16, 32

echo "scatter_count,count,start,end" > hpexome_timings.csv
for scatter_count in 1 2 4 8 16 32
do
    for count in 1 2 3 4 5
    do
        start="$(date +%s)"
        hpexome \
            --bam NA12878.sorted.rgfix.bam \
            --genome ucsc.hg19.fasta \
            --dbsnp dbsnp_138.hg19.vcf \
            --indels Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
            --indels 1000G_phase1.indels.hg19.sites.vcf \
            --sites 1000G_phase1.snps.high_confidence.hg19.sites.vcf \
            --sites 1000G_omni2.5.hg19.sites.vcf \
            --scatter_count ${scatter_count} \
            --job_runner GridEngine \
            result_files_${scatter_count}_${count}
        end="$(date +%s)"
        echo "${scatter_count},${count},${start},${end}" >> hpexome_timings.csv
    done
done
