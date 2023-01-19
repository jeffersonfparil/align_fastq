/////////////////////////////////////////////////////////////////////
// VARIANT CALLING //
/////////////////////////////////////////////////////////////////////

process INDEX_DEDUP {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir_reads
        val reference_genome
        val groupings
    output:
        val dir_reads
        val reference_genome
        val groupings
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}
    echo 'Index deduplicated alignments'
    parallel -j !{task.cpus} \
        picard BuildBamIndex \
            I={} \
        ::: $(ls *_deduped.bam)

    echo "Output:"
    echo "  (1/1) {group}_deduped.bai"
    '''
}

process ADD_READ_GROUPS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir_reads
        val reference_genome
        val groupings
    output:
        val dir_reads
        val reference_genome
        val groupings
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}
    echo 'Add read groups to the deduplicated alignments so that GATK4 works'
    echo 'And automatically index them too!'
    parallel -j !{task.cpus} \
        !{projectDir}/../scripts/add_read_groups_for_GATK4.sh \
            {1} \
        ::: $(ls *_deduped.bam)
    
    echo "Output:"
    echo "  (1/2) {group}_deduped_RG.bam"
    echo "  (2/2) {group}_deduped_RG.bai"
    '''
}

process CALL_VARIANTS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir_reads
        val reference_genome
        val groupings
    output:
        val dir_reads
        val groupings
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}
    CPU=!{task.cpus}
    MEM=$(echo !{task.memory} | cut -d' ' -f1)
    REF=!{reference_genome}
    
    echo "Call variants and automatically index the vcfs too!"
    for f in $(ls *_deduped_RG.bam)
    do
        name=$(echo $f | rev | cut -d'.' -f2- | rev)
        gatk --java-options "-Xmx${MEM}g" \
            HaplotypeCaller --native-pair-hmm-threads ${CPU} \
            -R ${REF} \
            -I ${f} \
            -O ${name}.vcf.gz \
            -ERC GVCF
    done

    echo "Output:"
    echo "  (1/2) {group}_deduped_RG.vcf.gz"
    echo "  (2/2) {group}_deduped_RG.vcf.gz.tbi"
    '''
}

process MERGE_VCFS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir_reads
        val groupings
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}
    CPU=!(task.cpus}
    MEM=$(echo !{task.memory} | cut -d' ' -f1)
    
    echo 'Find the groups of vcfs to merge.'
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

    echo "Merge vcfs"
    for g in $(cat group_names.tmp)
    do
        bcftools merge \
            --threads ${CPU} \
            -l vcflist-${g}.txt \
            -m all \
            -O z \
            -o ${g}.vcf.gz
        # gatk --java-options "-Xmx${MEM}g" \
        #     MergeVcfs \
        #     -I vcflist-${g}.txt \
        #     -O ${g}.vcf.gz 
    done

    echo "Cleanup"
    rm *.tmp

    echo "Output:"
    echo "  (1/2) {group}.vcf.gz"
    echo "  (2/2) vcflist-{group}.txt"
    '''
}

workflow {
    // INDEX_DEDUP(params.dir_reads, params.reference_genome, params.groupings) | \
    //     ADD_READ_GROUPS | \
    //     CALL_VARIANTS | \
        MERGE_VCFS(params.dir_reads, params.groupings)
}
