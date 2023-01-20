///////////////////////////////
// FILTER POOLS AND VARIANTS //
///////////////////////////////

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

    echo 'Filter synchronised pileups'
    for f in $(ls *.pileup)
    do
    done

    echo "Output:"
    echo "  (1/1) *_FILTERED.sync"
    '''
}

workflow {
    SYNCHRONISE(params.dir_reads,
                params.min_base_quality_Q)
}
