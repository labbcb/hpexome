# An automated workflow for processing whole-exome sequencing data

```bash
pip install --user --force git+https://github.com/labbcb/hpexome
```

Requirements

- BAM files must be sorted in `coordinate` mode. See [sort bam files](docs/sort_bam_files.sh) script.
- BAM files must have `@RG` tags with `ID, SM, LB, PL and PU` information. Seel [fix rg tag](docs/fix_rg_tag_bam.sh) script.

```bash
hpexome \
    --bam aln/NA12877_S1.sorted.rgfix.bam \
    --bam aln/NA12878_S1.sorted.rgfix.bam \
    --genome ref/ucsc.hg19.fasta  \
    --dbsnp ref/dbsnp_138.hg19.vcf \
    --indels ref/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    --indels ref/1000G_phase1.indels.hg19.sites.vcf \
    --sites ref/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
    --sites ref/1000G_omni2.5.hg19.sites.vcf \
    --unified_vcf \
    --num_data_threads 4 \
    --num_threads_per_data_thread 4 \
    res
```