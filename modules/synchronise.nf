///////////////////////////////////////////
// SYNCHRONISE THE PILED-UP MAPPED READS //
///////////////////////////////////////////

process SYNCHRONISE {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir_reads
        val min_base_quality_Q
    output:
        val dir_reads
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}
    MEM=$(echo !{task.memory} | cut -d' ' -f1)
    CPU=!{task.cpus}
    BAQ=!{min_base_quality_Q}

    echo 'Synchronise pileups'
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

    echo "Output:"
    echo "  (1/1) *.sync"
    '''
}

workflow {
    SYNCHRONISE(params.dir_reads,
                params.min_base_quality_Q)
}
