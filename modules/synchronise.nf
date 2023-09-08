///////////////////////////////////////////
// SYNCHRONISE THE PILED-UP MAPPED READS //
///////////////////////////////////////////

process SYNCHRONISE {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir_reads
        val groupings
        val min_base_quality_Q
        val popoolation2_not_poolgen
        val poolgen_min_coverage
        val poolgen_min_allele_freq
    output:
        val dir_reads
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}
    MEM=$(echo !{task.memory} | cut -d' ' -f1)
    CPU=!{task.cpus}
    BAQ=!{min_base_quality_Q}
    POPOOLATION_NOT_POOLGEN=!{popoolation2_not_poolgen}
    POOLGEN_MIN_COV=!{poolgen_min_coverage}
    POOLGEN_MIN_FRQ=!{poolgen_min_allele_freq}

    echo 'Synchronise pileups'
    if [ ${POPOOLATION_NOT_POOLGEN} == "TRUE" ] || [ ${POPOOLATION_NOT_POOLGEN} == "True" ] || [ ${POPOOLATION_NOT_POOLGEN} == "true" ] || [ ${POPOOLATION_NOT_POOLGEN} == "T" ] || [ ${POPOOLATION_NOT_POOLGEN} == "t" ]
    then
        for f in $(ls *.pileup)
        do
            out=$(echo $f | rev | cut -d"." -f2- | rev).sync
            java -ea -Xmx${MEM}g -jar \
            !{projectDir}/../scripts/popoolation2_1201/mpileup2sync.jar \
                --input ${f} \
                --output ${out} \
                --fastq-type sanger \
                --min-qual ${BAQ} \
                --threads ${CPU}
        done
    else
        cut -d',' -f1 !{groupings} | sort | uniq > group_names.tmp
        for g in $(cat group_names.tmp)
        do
            echo "#pop,pool_size,pheno" > mock_pheno-${g}.txt
            for l in $(cat bamlist-${g}.txt)
            do
                echo $(echo $(basename $l) | \
                    sed 's/.bam//g' | \
                    sed 's/paired-//g'),1,0 >> mock_pheno-${g}.txt
            done
            poolgen pileup2sync \
                -f ${g}.pileup \
                -p mock_pheno-${g}.txt \
                --phen-delim , \
                --phen-name-col 0 \
                --phen-pool-size-col 1 \
                --phen-value-col 2 \
                --min-coverage ${POOLGEN_MIN_COV} \
                --min-allele-frequency ${POOLGEN_MIN_FRQ} \
                --n-threads ${CPU} \
                -o ${g}.sync
        done
    fi

    echo "Output:"
    echo "  (1/1) *.sync"
    '''
}

workflow {
    SYNCHRONISE(params.dir_reads,
                params.groupings,
                params.min_base_quality_Q,
                params.popoolation2_not_poolgen,
                params.poolgen_min_coverage,
                params.poolgen_min_allele_freq)
}
