# An automated tool for processing whole-exome sequencing data

Whole-exome sequencing has been widely used in clinical applications for the identification of the genetic causes of several diseases.
_HPexome_ automates many data processing tasks for exome-sequencing data analysis of large-scale cohorts.
Given ready-analysis alignment files it is capable of breaking input data into small genomic regions to efficiently process in parallel on cluster-computing environments.
It relies on Queue workflow execution engine and GATK variant calling tool and its best practices to output high-confident unified variant calling file.
Our workflow is shipped as Python command line tool making it easy to install and use. 

Requirements

- BAM files must be sorted in `coordinate` mode. See [sort bam files](https://github.com/labbcb/hpexome/blob/master/docs/sort_bam_files.sh) script.
- BAM files must have `@RG` tags with `ID, SM, LB, PL and PU` information. See [fix rg tag](https://github.com/labbcb/hpexome/blob/master/docs/fix_rg_tag_bam.sh) script.

## Example

The following command line takes a list of ready-analysis BAM files stored in `alignment_files` directory and reference genomes files (version b37).
Then it breaks input data into smaller parts (`--scatter_count 16`) and submits to SGE batch system (`--job_runner PbsEngine`).
All samples will be merged into a single VCF files (`--unified_vcf`) and output files will be written in `result_files` directory.

```bash
hpexome \
    --bam alignment_files \
    --genome references/b37/human_g1k_v37_decoy.fasta  \
    --dbsnp references/b37/dbsnp_138.b37.vcf \
    --indels references/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
    --indels references/b37/1000G_phase1.indels.b37.vcf \
    --sites references/b37/1000G_phase1.snps.high_confidence.b37.vcf \
    --sites references/b37/1000G_omni2.5.b37.vcf \
    --unified_vcf \
    --scatter_count 16 \
    --job_runner GridEngine \
    result_fies
```
