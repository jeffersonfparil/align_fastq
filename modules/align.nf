/////////////////////////////////////////
// ALIGN READS TO THE REFERENCE GENOME //
/////////////////////////////////////////

process ALIGN {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val reference_genome
        val dir_reads
        val ext_read_1
        val ext_read_2
        val min_mapping_quality_Q
    output:
        val dir_reads
        val ext_read_1
        val ext_read_2
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}
    ext1=$(echo !{ext_read_1})
    ext2=$(echo !{ext_read_2})
    reference_genome_no_ext=$(echo !{reference_genome} | rev | cut -d'.' -f2- | rev)
    for r1 in $(ls paired-*${ext1})
    do
        r2=${r1%${ext1}*}${ext2}
        bamout=${r1%${ext1}*}.bam
        bwa mem \
            -t !{task.cpus} \
            ${reference_genome_no_ext} \
            ${r1} \
            ${r2} | \
        samtools view \
            -@ !{task.cpus} \
            -bu \
            -q !{min_mapping_quality_Q} \
            -T !{reference_genome} | \
        samtools sort \
            -@ !{task.cpus} \
            -o ${bamout} \
        || continue
    done
    
    echo "Output:"
    echo "  (1/1) paired-*.bam"
    '''
}

workflow {
    ALIGN(params.reference_genome,
          params.dir_reads,
          params.ext_read_1,
          params.ext_read_2,
          params.min_mapping_quality_Q)
}
