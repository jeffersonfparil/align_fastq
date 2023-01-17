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

## Per module (Nextflow compatabile with gatk3 is v21 and so needs -dsl2 flag to reconise params.* syntax)

```shell
time nextflow run -dsl2 modules/setup.nf -c config/params.config
time nextflow run -dsl2 modules/trim_and_qc.nf -c config/params.config
time nextflow run -dsl2 modules/align.nf -c config/params.config
time nextflow run -dsl2 modules/pileup.nf -c config/params.config

time nextflow run -dsl2 modules/synchronise.nf -c config/params.config

time nextflow run -dsl2 modules/dedup.nf -c config/params.config
time nextflow run -dsl2 modules/variant_calling.nf -c config/params.config



```