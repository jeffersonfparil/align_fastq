/////////////////////////////
// PILEUP THE MAPPED READS //
/////////////////////////////

process PILEUP {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val reference_genome
        val dir_reads
        val groupings
        val min_mapping_quality_Q
        val min_base_quality_Q
        val window_size
        val n_chromosomes
    output:
        val dir_reads
        val window_size
        val n_chromosomes
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}

    echo 'Find the groups of alignments (*.bam) to pileup.'
    cut -d',' -f1 !{groupings} | sort | uniq > group_names.tmp
    for g in $(cat group_names.tmp)
    do
        grep "^${g}," !{groupings} | cut -d',' -f2 | sort | uniq > names_within_group.tmp
        rm bamlist-${g}.txt || touch bamlist-${g}.txt
        for name in $(cat names_within_group.tmp)
        do
            find !{dir_reads} -name '*.bam' | grep -v "unpaired-" | grep -v "deduped" | grep "${name}" | head -n1 >> bamlist-${g}.txt
        done
    done
    
    echo 'Pileup in parallel.'
    parallel \
        -j !{task.cpus} \
        --tmpdir !{dir_reads} \
        sh !{projectDir}/../scripts/pileup.sh \
            {} \
            !{reference_genome} \
            !{min_mapping_quality_Q} \
            !{min_base_quality_Q} \
        ::: $(ls bamlist-*.txt)

    echo "Clean-up"
    rm *.tmp

    echo "Output:"
    echo "  (1/2) bamlist-*.txt"
    echo "  (2/2) *.pileup"
    '''
}

process COVERAGE_BREADTH_AND_DEPTH {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir_reads
        val window_size
        val n_chromosomes
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}

    ### Estimate summary statistics for the breadth and depth of coverage
    parallel -j !{task.cpus} \
        julia !{projectDir}/../scripts/pileup_stats.jl \
            {} \
            !{window_size} \
            !{n_chromosomes} \
    ::: $(ls *.pileup)

    echo "Output:"
    echo "  (1/1) *-coverage_stats.png"
    '''
}

workflow {
    PILEUP(params.reference_genome,
           params.dir_reads,
           params.groupings,
           params.min_mapping_quality_Q,
           params.min_base_quality_Q,
           params.window_size,
           params.n_chromosomes) | \
    COVERAGE_BREADTH_AND_DEPTH
}
