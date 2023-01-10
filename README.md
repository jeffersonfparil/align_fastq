# align_fastq
Alignment or mapping of nucleotide sequencing reads

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
conda activate align_fastq
```

## Quickstart

Run the example pipeline:
```shell
cd align_fastq/
nano config/params.config # replace dir = '/data-weedomics-3/TEST_PSEUDOMONAS' with a valid path in your machine
chmod +x run.sh
time ./run.sh
```

## Per module

```shell
time nextflow run modules/trim_and_qc.nf -c config/params.config
time nextflow run modules/setup_reference_genome.nf -c config/params.config
time nextflow run modules/align.nf -c config/params.config
time nextflow run modules/pileup.nf -c config/params.config
```