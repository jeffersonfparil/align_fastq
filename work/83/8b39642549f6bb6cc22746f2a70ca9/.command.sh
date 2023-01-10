#!/usr/bin/env bash
echo 'Fix the reference genome format.'
mv /data-weedomics-1/align_fastq/test/ref/Lolium_rigidum_genome.fasta /data-weedomics-1/align_fastq/test/ref/Lolium_rigidum_genome.fasta.bk
python3 /data-weedomics-1/align_fastq/modules/../scripts/fix_Reference_genome_format.py         /data-weedomics-1/align_fastq/test/ref/Lolium_rigidum_genome.fasta.bk         /data-weedomics-1/align_fastq/test/ref/Lolium_rigidum_genome.fasta
