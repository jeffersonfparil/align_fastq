////////////////////////////////////////////////////////////////////
// REMOVE DUPLICATE READS FOR INDIVIDUAL GENOTYPE VARIANT CALLING //
////////////////////////////////////////////////////////////////////

process DEDUP {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    echo "Deduplication with Picard"
    parallel -j !{task.cpus} \
        !{projectDir}/../scripts/dedup.sh \
            {} \
        ::: $(ls *.bam)

    echo "Output:"
    echo "  (1/2) *_deduped.bam"
    echo "  (2/2) *_dedup_metrics.txt"
    '''
}

workflow {
    DEDUP(params.dir_reads)
}
