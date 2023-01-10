#!/usr/bin/env bash
cd /data-weedomics-1/align_fastq/test/reads
ext1=$(echo _R1.fastq.gz)
ext2=$(echo _R2.fastq.gz)

for r1 in $(ls *${ext1})
do
    r2=${r1%${ext1}*}${ext2}
    trimmomatic             PE             -threads 32             ${r1} ${r2}             paired-${r1} paired-${r2}             unpaired-${r1} unpaired-${r2}             ILLUMINACLIP:/data-weedomics-1/align_fastq/test/IDT_for_Illumina_TruSeq_UD_and_CD_indexes.fa:2:30:10:2:True             LEADING:3             TRAILING:3             MINLEN:36 || continue
done


echo "Output:"
echo "  (1/2) {ORTHONAME}.fasta"
echo "  (2/2) {ORTHONAME}-{SPECIES}.fasta"
