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

    for r1 in $(ls *${ext1})
    do
        r2=${r1%${ext1}*}${ext2}
        trimmomatic \
            PE \
            -threads !{task.cpus} \
            ${r1} ${r2} \
            paired-${r1} unpaired-${r1} \
            paired-${r2} unpaired-${r2} \
            ILLUMINACLIP:!{adapters}:2:30:10:2:True \
            LEADING:3 \
            TRAILING:3 \
            MINLEN:36 \
        || continue
    done
    
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
