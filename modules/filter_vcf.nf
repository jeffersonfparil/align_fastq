/////////////////////////////////
// FILTER SAMPLES AND VARIANTS //
/////////////////////////////////

process FILTER {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir_reads
        val groupings
        val call_quality
        val minimum_depth
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}
    CPU=!(task.cpus}
    CALLQ=!{min_call_quality}
    DEPTH=!{min_depth}

    # echo "NOTE: QUAL = -10*log10(probability of homozygous reference genotype), e.g. QUAL=10 means there is 10% chance that there are no variants in the site, while QUAL=20 means 1% and QUAL=30 means 0.1%."
    echo "NOTE: GQ = -10*log10(probability of incorrect variant call), e.g. QUAL=10 means there is 10% chance that call was incorrect, while QUAL=20 means 1% and QUAL=30 means 0.1%."
    echo "NOTE: DP is the depth of sequencing coverage."

    echo 'Find the names of the groups of merged vcfs.'
    cut -d',' -f1 !{groupings} | sort | uniq > group_names.tmp
    for g in $(cat group_names.tmp)
    do
        grep "^${g}," !{groupings} | cut -d',' -f2 | sort | uniq > names_within_group.tmp
        rm vcflist-${g}.txt || touch vcflist-${g}.txt
        for name in $(cat names_within_group.tmp)
        do
            find !{dir_reads} -name '*_deduped_RG.vcf.gz' | grep "paired-" | grep "${name}" | head -n1 >> vcflist-${g}.txt
        done
    done

    echo "Filter by minimum variant calling error (GQ) and minimum sequencing depth (DP)."
    for g in $(cat group_names.tmp)
    do
        vcf=${g}.vcf.gz
        bcftools view \
            --threads ${CPU} \
            -e "GQ >= ${CALLQ} || DP > ${DEPTH}" \
            ${vcf} \
            -O z \
            -o ${g}-FILTERED_GQ${CALLQ}_DP${DEPTH}.vcf.gz
    done

    echo "Output:"
    echo "  (1/1) {group}-FILTERED_GQ${CALLQ}_DP${DEPTH}.vcf.gz"
    '''
}

workflow {
    FILTER(params.dir_reads, params.groupings, params.min_call_quality, params.min_depth) | \
}
