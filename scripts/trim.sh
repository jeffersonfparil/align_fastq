#!/usr/bin/env bash
r1=$1
ext1=$2
ext2=$3
adapters=$4
r2=${r1%${ext1}*}${ext2}

trimmomatic \
    PE \
    -threads 1 \
    ${r1} ${r2} \
    paired-${r1} unpaired-${r1} \
    paired-${r2} unpaired-${r2} \
    ILLUMINACLIP:${adapters}:2:30:10:2:True \
    LEADING:3 \
    TRAILING:3 \
    MINLEN:36 \
    || exit