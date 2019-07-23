# DOC: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_AddOrReplaceReadGroups.php

# Path to Picard tool
PICARD=picard.jar

for file in "$@"
do
    filename="$(basename $file)"
    name=${filename%_S1*}
    java -jar $PICARD AddOrReplaceReadGroups \
        I=$file \
        O="${file%.*}".rgfix.bam \
        RGID=$name \
        RGSM=$name \
        RGLB=platinum \
        RGPL=illumina \
        PU=unit1 \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT # this tells Picard to ignore "MAPQ should be 0 for unmapped read"
done