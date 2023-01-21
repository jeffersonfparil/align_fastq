# align_fastq
Workflow to align nucleotide sequencing data (i.e. short (<1kb) paired-end reads) into reference genomes, call and count variants.

## Installation

1. Dowload [align_fastq](https://github.com/jeffersonfparil/align_fastq.git) repository
```shell
git clone https://github.com/jeffersonfparil/align_fastq.git
```

2. Install conda
```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh ./Miniconda3-latest-Linux-x86_64.sh
```

3. **Log out of the current session, then log back in** to finish conda setup.

4. Import and activate [align_fastq conda environment](align_fastq.yml)
```shell
conda env create -n align_fastq --file align_fastq/align_fastq.yml
conda env create -n align_fastq --file align_fastq/align_fastq_bcftools.yml
conda activate align_fastq
```

## Configuration

Setup your alignment, and variant calling/counting pipeline. You will find the following files in [align_fastq/config](config/) folder.

1. [`groupings.txt`](config/groupings.txt): list of intended groupings of the paired-end reads, i.e. tells how to merge reads into vcf, pileup, and/or sync files

    - Formatted as headerless, two-columned, comma-delimited file
    - Column 1: group name, i.e. the name of the final `*.vcf.`gz and/or `*.sync` file/s
    - Column 2: base name of the paired-end read excluding the read 1 or 2 identifier (see item 2 above for more details) which belongs to the correspponding group in column 1
    - Note: one or more paired-end reads can belong to a single group

2. [`params.config`](config/params.config): list of parameters specific to your data set and variant calling (i.e. for individual sequencing data) or variant counting (i.e. for pool sequencing data)

    - **dir_reads**: location of the short paired-end sequencing reads (e.g. '/data-weedomics-1/align_fastq/test/reads')
    - **ext_read_1**: suffix of the read 1 (e.g. '_R1.fastq.gz'; note that reads of each paired-end read pair needs to have the same base name excluding the read 1 or 2 identifier)
    - **ext_read_2**: suffix of the read 2 (e.g. '_R2.fastq.gz')
    - **adapters**: fasta file containing adapter sequences used during sequencing (e.g. "${projectDir}/../test/IDT_for_Illumina_TruSeq_UD_and_CD_indexes.fa")
    - **reference_genome**: reference genome in fasta format (e.g. "${projectDir}/../test/ref/Lolium_rigidum_genome.fasta")
    - **min_mapping_quality_Q**: minimum mapping quality (e.g. 20 which equates to 0.01 error rate)
    - **min_base_quality_Q**: minimum sequenced base quality (e.g. 20 which equates to 0.01 error rate)
    - **groupings**: list of intended groupings of the paired-end reads (e.g. "${projectDir}/../config/groupings.txt")

3. [`process.config`](config/process.config): list of the computing resource allocation availble to you. Assign the number of cpus and memory capacity to use for low and high resources tasks:

    - LOW_MEM_LOW_CPU
    - HIGH_MEM_HIGH_CPU

## Run

Run for individual genotyping sequence data:
```shell
cd align_fastq/
time ./run_indiseq.sh
```

Run for pool genotyping sequence data:
```shell
cd align_fastq/
time ./run_poolseq.sh
```

Run each module individually, usually to troubleshoot:
```shell
time nextflow run modules/setup.nf -c config/params.config              ### Setup reference genome and Julia packages
time nextflow run modules/trim_and_qc.nf -c config/params.config        ### Remove adapters and perform quality check of the raw reads
time nextflow run modules/align.nf -c config/params.config              ### Align the reads to the reference genome
time nextflow run modules/pileup.nf -c config/params.config             ### Pileup the reference genome and assess the distribution of the breadth and depth of sequencing
### For Pool-seq
time nextflow run modules/synchronise.nf -c config/params.config        ### For Pool-seq data: convert the pileups into syncrhonised pileup format
### For Indi-seq
time nextflow run modules/dedup.nf -c config/params.config              ### For Indi-seq data: remove PCR duplicates from the raw reads
time nextflow run modules/variant_calling.nf -c config/params.config    ### For Indi-seq data: perform varant calling but first index the alignments, and add read groups, and finally merge the VCF files
```
