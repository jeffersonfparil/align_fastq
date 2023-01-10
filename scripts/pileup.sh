#!/usr/bin/env bash
BAMLIST=$1
REF=$2
MAPQ=$3
BASQ=$4
PILEUPOUT=$(echo ${BAMLIST} | sed 's/bamlist-//g' | sed 's/.txt//g').pileup

samtools mpileup \
    --min-MQ ${MAPQ} \
    --min-BQ ${BASQ} \
    --fasta-ref ${REF} \
    --bam-list ${BAMLIST} \
    --output ${PILEUPOUT}
