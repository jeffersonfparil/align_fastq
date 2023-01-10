#!/usr/bin/env bash
cd /data-weedomics-1/align_fastq/test/reads
ext1=$(echo _R1.fastq.gz)
ext2=$(echo _R2.fastq.gz)

for f in $(ls *${ext1}) $(ls *${ext2})
do
    fqc             -q ${f} > fastqc-${f}.html         || continue
done

echo "Output:"
echo "  (1/1) fastqc-*.html"
