#!/usr/bin/env bash
echo 'Index the reference genome for bwa.'
reference_genome_no_ext=$(echo /data-weedomics-1/align_fastq/test/ref/Lolium_rigidum_genome.fasta | rev | cut -d'.' -f2- | rev)
bwa index         -p ${reference_genome_no_ext}         -a bwtsw         /data-weedomics-1/align_fastq/test/ref/Lolium_rigidum_genome.fasta
