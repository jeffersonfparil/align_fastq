#!/bin/bash

echo "###########################################"
echo "ALIGN AND ALLELE COUNTING FOR POOL-SEQ DATA"
echo "###########################################"

echo "Setup reference genome and Julia packages"
time nextflow run modules/setup.nf -c config/params.config

echo "Remove adapters and perform quality check of the raw reads"
time nextflow run modules/trim_and_qc.nf -c config/params.config

echo "Align the reads to the reference genome"
time nextflow run modules/align.nf -c config/params.config

echo "Pileup the reference genome and assess the distribution of the breadth and depth of sequencing"
time nextflow run modules/pileup.nf -c config/params.config

echo "Syncrhonised the pileup files"
time nextflow run modules/synchronise.nf -c config/params.config
