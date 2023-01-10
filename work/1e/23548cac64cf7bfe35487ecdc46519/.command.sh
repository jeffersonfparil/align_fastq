#!/usr/bin/env bash
cd /data-weedomics-1/align_fastq/test/reads

echo 'Find the groups of alignments (*.bam) to pileup.'
cut -d',' -f1 /data-weedomics-1/align_fastq/modules/../config/groupings.txt | sort | uniq > group_names.tmp
for g in $(cat group_names.tmp)
do
    grep "^${g}," /data-weedomics-1/align_fastq/modules/../config/groupings.txt | cut -d',' -f2 | sort | uniq > names_within_group.tmp
    touch bamlist-${g}.txt
    for name in $(cat names_within_group.tmp)
    do
        ls | grep "paired-" | grep "${name}" | head -n1 >> bamlist-${g}.txt
    done
done

echo 'Pileup in parallel.'
parallel         -j 32         /data-weedomics-1/align_fastq/modules/../scripts/pileup.sh             {}             /data-weedomics-1/align_fastq/test/ref/Lolium_rigidum_genome.fasta             20             20         ::: $(ls bamlist-*.txt)

echo "Output:"
echo "  (1/2) bamlist-*.txt"
echo "  (2/2) *.pileup"
