#!/usr/bin/env bash
cd /data-weedomics-1/align_fastq/test/reads
ext1=$(echo _R1.fastq.gz)
ext2=$(echo _R2.fastq.gz)
reference_genome_no_ext=$(echo /data-weedomics-1/align_fastq/test/ref/Lolium_rigidum_genome.fasta | rev | cut -d'.' -f2- | rev)
for r1 in $(ls paired-*${ext1})
do
    r2=${r1%${ext1}*}${ext2}
    bamout=${r1%${ext1}*}.bam
    bwa mem             -t 32             ${reference_genome_no_ext}             ${r1}             ${r2} |         samtools view             -@ 32             -b             -q 20             -T /data-weedomics-1/align_fastq/test/ref/Lolium_rigidum_genome.fasta |         samtools sort             -@ 32             -o ${bamout}         || continue
done

echo "Output:"
echo "  (1/1) paired-*.bam"
