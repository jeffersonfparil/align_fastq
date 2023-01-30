# Installation from scratch

1. Install conda
```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh ./Miniconda3-latest-Linux-x86_64.sh
echo "To finish the conda setup - log-out then log back in."
```
2. Configure conda
```shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. Create a new conda environment
```shell
conda create --name "align_fastq"
conda env list
conda activate align_fastq
```

4. Install software
```shell
conda install -y parallel
conda install -y -c bioconda nextflow
conda install -y -c bioconda trimmomatic
conda install -y -c bioconda fastqc-rs
conda install -y -c bioconda biopython
conda install -y -c bioconda bwa
conda install -y -c bioconda 'samtools>=1.10'
conda install -y -c bioconda java-jdk
conda install -y -c bioconda picard
conda install -y -c bioconda gatk4
conda install -y -c bioconda bcftools
conda install -y -c bioconda 'julia=1.7.1'
```

5. Export environment and create a new environment based on the exported settings
```shell
conda env export -n align_fastq > align_fastq.yml
```

6. Import conda environment
```shell
conda env create -n align_fastq_import --file align_fastq.yml
```

7. Dowload this comparative_genomics repository
```shell
git clone https://github.com/jeffersonfparil/align_fastq.git
```