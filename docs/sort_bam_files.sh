# USAGE bash sort_bam_files.sh file.bam [file2.bam ...]
# Doc: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_SortSam.php

# Path to Picard tool
PICARD=picard.jar

for file in "$@"
do
    java -jar $PICARD SortSam \
        INPUT=$file \
        OUTPUT="${file%.*}".sorted.bam \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT # this tells Picard to ignore "MAPQ should be 0 for unmapped read"
done
