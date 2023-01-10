#!/bin/bash
DIR=$1
cd $DIR
echo "TEST" > test-2.tmp
date >> test-2.tmp


# Andhika's Variant calling protocol

# 1. Align the mother plants to the new genome


for SRID in `ls /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/AU**.gz | xargs -n 1 basename | cut -d_ -f1-2`
do bwa mem ~/Clim_GWAS/RAGWEED/NEW_GENOME/mon2034_05-30-2021_hic_output_fixed_header.fasta /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/${SRID}_1.fq.gz /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/${SRID}_2.fq.gz > /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/${SRID}_pe.newgen.sam -R $(echo '@RG\tID:'"${SRID}"'\tSM:bar')
done 



# 2. Convert SAM to BAM and sort by scaffold name

for SRID in `ls /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/*.sam | xargs -n 1 basename | cut -d_ -f1-2`
do samtools view -b ~/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/${SRID}_pe.newgen.sam > ~/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/$SRID.newgen.bam
done


for SRID in `ls /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/*.sam | xargs -n 1 basename | cut -d_ -f1-2`
do samtools sort ${SRID}.newgen.bam -o ${SRID}_sort.newgen.bam
done

# 3. Remove duplicates

ls /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/*.sam | xargs -n 1 basename | cut -d_ -f1-2 | parallel java -jar ~/Clim_GWAS/jars/picard.jar MarkDuplicates --REMOVE_DUPLICATES true -I "~/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/{1}_sort.newgen.bam" "-O ~/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/{1}_markdup.newgen.bam" "-M ~/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/{1}_metrics.newgen.txt"


# 4. index the bam files

cat /home/student.unimelb.edu.au/andhikap/Clim_GWAS/sratoolkit.2.10.9-ubuntu64/bin/x* | parallel java -jar ~/Clim_GWAS/picard.jar BuildBamIndex \
      I="~/Clim_GWAS/RAGWEED/SEQDATA/{1}_markdup.newgen.bam"

# 5. Genotyping using HaplotypeCaller in GATK4

for SRID in `ls /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/*.sam | xargs -n 1 basename | cut -d_ -f1-2`
do ~/Clim_GWAS/gatk-4.2.1.0/gatk HaplotypeCaller -R ~/Clim_GWAS/RAGWEED/ragweed2021/ragweed2021.fasta -I ~/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/${SRID}_markdup.newgen.bam -O ~/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/${SRID}.newgen.g.vcf.gz -ERC GVCF
done


# 6. Making the sample map file

ls /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/*.sam | xargs -n 1 basename | cut -d_ -f1-2 > sample_names.list

ls /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/*.newgen.g.vcf.gz > sample_files.list

paste -d"\t" sample_names.list sample_files.list > cohort.sample_map


# 7. Import the genotyped samples to genomicsdbimport and perform joint genotyping
parallel --jobs 5 ./genomicsdbimport_parallel_script.sh $1 ::: {1..451}

	# for SRID in {1..451}
	# do

	# ~/Clim_GWAS/gatk-4.2.1.0/gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport \
	#        --genomicsdb-workspace-path /home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/GENOMICSDBIMPORT_SPACE_$1 \
	#        --batch-size 20 \
	#        --sample-name-map cohort.sample_map \
	#        --tmp-dir ~/Clim_GWAS/RAGWEED/GENOMICSDB_TEMPDIR \
	#        --reader-threads 5 \
	#        --L ScBFxKa_$1

	#        ~/Clim_GWAS/gatk-4.2.1.0/gatk GenotypeGVCFs --allow-old-rms-mapping-quality-annotation-data -R ~/Clim_GWAS/RAGWEED/ragweed2021/ragweed2021.fasta -V gendb:///home/student.unimelb.edu.au/andhikap/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/GENOMICSDBIMPORT_SPACE_$1 -O ~/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/NEWGEN_TEMP_VCFs/$1.L7.vcf.gz

	#         rm -rf ~/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/GENOMICSDBIMPORT_SPACE_$1
	# done

# 8. Merge genotyped samples

ls ~/Clim_GWAS/RAGWEED/AZENTA_FASTQ_251122/N2213286_30-784346616_Lib_2022-11-25/221122-E00516A_L007/NEWGEN_TEMP_VCFs/*.vcf.gz >> inputs.newgen.list; sort -V  inputs.newgen.list > temp.list ; mv temp.list inputs.newgen.list


~/Clim_GWAS/gatk-4.2.0.0/gatk MergeVcfs -I inputs.newgen.list  -O fieldwork_combined.vcf.gz 