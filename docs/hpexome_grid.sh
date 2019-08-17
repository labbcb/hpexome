#!/bin/bash

for sc in 1 2 4 8 16
do
    time_total=0
    for count in 1 2 3 4 5
    do
        outdir="outdir_1bam_${sc}sc_${count}"
        mkdir $outdir
        time_before="$(date +%s)"
        java -jar Queue.jar -S Script.scala -R genome.fa -I input/A01.bam.rmdup.bam -known ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf -known ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf -knownSites clinvar_20140604.vcf -dbsnp clinvar_20140604.vcf -stand_call_conf 40 -stand_emit_conf 10 -minPruning 3 -l DEBUG -run -scattercount $sc -outdir $outdir -jobRunner GridEngine -jobEnv "make 1" -jobQueue all.q@lmbs045.cnptia.embrapa.br -| qsub -V -cwd -N "try_1bam_${sc}sc_${count}"
        time_after="$(date +%s)"
        ((time_elapsed=time_after - time_before))
        ((time_total+=time_elapsed))
        echo "scattercount = $sc count = $count elapsed = $time_elapsed"
    done
    ((mean=time_total / 5))
    echo "scattercount = $sc mean elapsed = $mean"
    echo
done

