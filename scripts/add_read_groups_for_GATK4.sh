#!/usr/bin/env bash
BAM=$1
NAME=$(echo ${BAM} | rev | cut -d'.' -f2- | rev)

picard \
    AddOrReplaceReadGroups \
        I=${BAM} \
        O=${NAME}_RG.bam \
        SORT_ORDER=coordinate \
        RGID=${NAME} \
        RGLB=NEBNext_Illumina \
        RGPL=illumina \
        RGSM=${NAME} \
        RGPU=wahevah \
        CREATE_INDEX=True
