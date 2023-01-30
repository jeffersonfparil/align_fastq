#!/usr/bin/env bash
f=$1
fqc \
    -q ${f} \
    > fastqc-${f}.html \
|| continue