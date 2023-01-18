#!/bin/bash

echo "###########################################"
echo "ALIGN AND VARIANT CALLING FOR INDI-SEQ DATA"
echo "###########################################"

echo "Setup reference genome and Julia packages"
time nextflow run modules/setup.nf -c config/params.config

echo "Remove adapters and perform quality check of the raw reads"
time nextflow run modules/trim_and_qc.nf -c config/params.config

echo "Align the reads to the reference genome"
time nextflow run modules/align.nf -c config/params.config

echo "Pileup the reference genome and assess the distribution of the breadth and depth of sequencing"
time nextflow run modules/pileup.nf -c config/params.config

echo "Remove PCR duplicates from the raw reads"
time nextflow run modules/dedup.nf -c config/params.config

echo "Perform varant calling but first index the alignments, and add read groups, and finally merge the VCF files"
time nextflow run modules/variant_calling.nf -c config/params.config
