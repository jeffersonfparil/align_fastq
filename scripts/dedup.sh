#!/usr/bin/env bash
f=$1
name=$(echo $f | rev | cut -d"." -f2- | rev)
picard MarkDuplicates \
    --REMOVE_DUPLICATES true \
    -I ${f} \
    -O ${name}_deduped.bam \
    -M ${name}_dedup_metrics.txt