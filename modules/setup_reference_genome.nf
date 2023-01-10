////////////////////////////////
// SETUP THE REFERENCE GENOME //
////////////////////////////////

process FIX {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val reference_genome
    output:
        val reference_genome
    shell:
    '''
    #!/usr/bin/env bash
    
    echo 'Fix the reference genome format.'
    mv !{reference_genome} !{reference_genome}.bk
    python3 !{projectDir}/../scripts/fix_reference_genome_format.py \
        !{reference_genome}.bk \
        !{reference_genome}
    
    echo "Output:"
    echo "  (1/2) {reference_genome}"
    echo "  (2/2) {reference_genome}.bk"
    '''
}

process BWA_INDEX {
    label "LOW_MEM_LOW_CPU"
    input:
        val reference_genome
    shell:
    '''
    #!/usr/bin/env bash
    
    echo 'Index the reference genome for bwa.'
    reference_genome_no_ext=$(echo !{reference_genome} | rev | cut -d'.' -f2- | rev)
    bwa index \
        -p ${reference_genome_no_ext} \
        -a bwtsw \
        !{reference_genome}
    
    echo "Output:"
    echo "  (1/3) {reference_genome_no_ext}.amb"
    echo "  (2/3) {reference_genome_no_ext}.ann"
    echo "  (3/3) {reference_genome_no_ext}.pac"
    '''
}

process SAMTOOLS_INDEX {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val reference_genome
    shell:
    '''
    #!/usr/bin/env bash
    
    echo 'Index the reference genome for samtools.'
    samtools faidx \
        !{reference_genome} \
        -@ !{task.cpus}
   
    echo "Output:"
    echo "  (1/1) {reference_genome_no_ext}.sa"
    '''
}

workflow {
    FIX(params.reference_genome) | \
    (BWA_INDEX & SAMTOOLS_INDEX)
}