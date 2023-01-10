#!/usr/bin/env bash

echo 'Index the reference genome for samtools.'
samtools faidx         /data-weedomics-1/align_fastq/test/ref/Lolium_rigidum_genome.fasta         -@ 32

echo "Output:"
echo "  (1/1) {reference_genome}.fai"
