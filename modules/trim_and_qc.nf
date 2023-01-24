///////////////////////////////////
// TRIM READS AND ASSESS QUALITY //
///////////////////////////////////

process TRIM {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir_reads
        val ext_read_1
        val ext_read_2
        val adapters
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
    adapters=!{adapters}

    parallel -j !{taks.cpus} \
        !{projectDir}/../scripts/trim.sh \
            {} \
            ${ext1} \
            ${ext2} \
            ${adapters} \
        ::: $(ls *${ext1})
    
    echo "Output:"
    echo "  (1/2) paired-*"
    echo "  (2/2) unpaired-*"
    '''
}

process QC {
    label "LOW_MEM_LOW_CPU"
    input:
        val dir_reads
        val ext_read_1
        val ext_read_2
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir_reads}
    ext1=$(echo !{ext_read_1})
    ext2=$(echo !{ext_read_2})

    for f in $(ls *${ext1}) $(ls *${ext2})
    do
        fqc \
            -q ${f} > fastqc-${f}.html \
        || continue
    done
    
    echo "Output:"
    echo "  (1/1) fastqc-*.html"
    '''
}

workflow {
    TRIM(params.dir_reads,
         params.ext_read_1,
         params.ext_read_2,
         params.adapters) | \
    QC
}
