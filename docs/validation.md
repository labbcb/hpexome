---
title: "Validating hpexome tool"
output: 
  html_document: 
    keep_md: yes
---

Command line executed on EMBRAPA grid:

```bash
time java -Djava.io.tmpdir=.hpexome/ -jar HPexome.jar -R hpexome_validation/genome.chr.fa \
-I hpexome_validation/NA12877_S1.sort.fix.bam \
-known hpexome_validation/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.chr.vcf \
-known hpexome_validation/ALL.wgs.low_coverage_vqsr.20101123.indels.sites.chr.vcf \
-knownSites hpexome_validation/clinvar_20140604.chr.vcf -dbsnp hpexome_validation/clinvar_20140604.chr.vcf \
-stand_call_conf 40 -stand_emit_conf 10 -minPruning 3 -l DEBUG \
-run -nct 8 -nt 8 -scattercount 8 \
-jobRunner GridEngine -jobEnv "make 1" \
-outdir hpexome_validation/validate_NA12877_try1 | qsub -V -cwd -N validate_NA12877_try1
```

Download VCF reference: ftp://ussd-ftp.illumina.com/hg19/8.0.1/NA12877/

Packages

```r
library(VariantAnnotation)
library(knitr)
library(ggplot2)
```


Normalize VCF files:

```r
workDir <- "/home/bioinf/hpexome_validation"
genomeFile <- file.path(workDir, "genome.chr.fa")
NA12877.ref <- file.path(workDir, "NA12877.vcf.gz")
NA12877.ref.vt <- file.path(workDir, "NA12877.vt.vcf")
NA12877.vcf <- file.path(workDir, "validate_NA12877_try1/NA12877_S1.sort.fix.HC.raw.vt.vcf.bz")
NA12877.vcf.vt <- file.path(workDir, "validate_NA12877_try1/NA12877_S1.sort.fix.HC.raw.vt.vcf")
```

    vt decompose /home/bioinf/hpexome_validation/NA12877.vcf.gz | vt normalize -r /home/bioinf/hpexome_validation/genome.chr.fa - | vt uniq - -o /home/bioinf/hpexome_validation/NA12877.vt.vcf
    vt decompose /home/bioinf/hpexome_validation/validate_NA12877_try1/NA12877_S1.sort.fix.HC.raw.vt.vcf.bz | vt normalize -r /home/bioinf/hpexome_validation/genome.chr.fa - | vt uniq - -o /home/bioinf/hpexome_validation/validate_NA12877_try1/NA12877_S1.sort.fix.HC.raw.vt.vcf


```r
NA12877.ref.vt.bz <- paste0(NA12877.ref.vt, ".bz")
# bgzip(NA12877.ref.vt, NA12877.ref.vt.bz)
NA12877.vcf.vt.bz <- paste0(NA12877.vcf.vt, ".bz")
# bgzip(NA12877.ref.vt, NA12877.vcf.vt.bz)
```

Load VCF files into R:

```r
genome <- "hg19"
ref <- readVcf(NA12877.ref.vt.bz, genome)
vcf <- readVcf(NA12877.vcf.vt.bz, genome)
```

Filter VCFs by selected chromosomes

```r
chr.select <- paste0("chr", 1:22)
seqlevels(ref, force=TRUE) <- seqlevels(ref)[seqlevels(ref) %in% chr.select]
seqlevels(vcf, force=TRUE) <- seqlevels(vcf)[seqlevels(vcf) %in% chr.select]
```

Remove variants from reference that have genotype "0/."

```r
idx <- which(!grepl("0(\\/|\\|)\\.", geno(ref)$GT))
ref <- ref[idx, ]
```

Number of variants:

```r
kable(data.frame(reference=nrow(ref), hpexome=nrow(vcf)))
```



| reference| hpexome|
|---------:|-------:|
|   4425657| 4706967|

Subset VCF files using common variants between the two VCFs:

```r
vcf.ranges <- rowRanges(vcf)
ref.ranges <- rowRanges(ref)
overlap <- findOverlaps(vcf.ranges, ref.ranges)
vcf.sub <- vcf[queryHits(overlap), ]
ref.sub <- ref[subjectHits(overlap), ]
```

Number of variants after subset (must be the same value for the two file):

```r
kable(data.frame(reference=nrow(ref.sub), hpexome=nrow(vcf.sub)))
```



| reference| hpexome|
|---------:|-------:|
|   4356821| 4356821|

The hpexome tool not found 68836 variants.
However this tool found 350146 new variants.
Because hpxome found 281310 more variants than reference.

Create genotype count table

```r
vcf.gt <- geno(vcf.sub)$GT
ref.gt <- sub("\\|", "/", geno(ref.sub)$GT)
ref.gt <- sub("1\\/0", "0/1", ref.gt)

genoCount <- table(ref.gt, vcf.gt)
genoCount
```

```
##       vcf.gt
## ref.gt       .     0/1     1/.     1/1
##    .     34160    7871   34178    8255
##    0/1    2485 2584861    3347    4205
##    1/.   30938    2958   30945    2539
##    1/1     642    4807     663 1603967
```

Degree of concordance

```r
sum(diag(genoCount)) / sum(genoCount) * 100
```

```
## [1] 97.63846
```

Create intevals for `QG` value:

```r
gq <- geno(vcf.sub)$GQ
breaks <- unique(quantile(gq, seq(0, 1, .01), na.rm=TRUE))
intervals <- cut(gq, breaks)
```

This function count the calculates the number of correct genotypes:

```r
f <- function(interval, vcf.gt, ref.gt, intervals)
{
  idx <- interval == intervals
  sum(vcf.gt[idx] == ref.gt[idx], na.rm=TRUE) / length(ref.gt[idx]) * 100
}
```

Create table for intervals and correct genotypes:

```r
correct <- sapply(levels(intervals), f, vcf.gt, ref.gt, intervals, USE.NAMES = FALSE)
df <- data.frame(intervals=levels(intervals), correct)
kable(df)
```



|intervals |  correct|
|:---------|--------:|
|(0,42]    | 18.59625|
|(42,64]   | 21.80658|
|(64,80]   | 22.67760|
|(80,92]   | 23.52834|
|(92,99]   | 96.33104|

Plot intervals

```r
ggplot(df, aes(x=intervals, y=correct)) +
    geom_bar(stat="identity") +
    labs(x="QG value intervals", y="% correct variant")
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png) 
