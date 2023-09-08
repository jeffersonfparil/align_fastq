///////////////////////////////////////////////////
// SETUP THE REFERENCE GENOME AND JULIA PACKAGES //
///////////////////////////////////////////////////

process GENOME_FIX {
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

process GENOME_BWA_INDEX {
    errorStrategy 'terminate'
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

process GENOME_SAMTOOLS_INDEX {
    errorStrategy 'terminate'
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

process GENOME_GATK4_INDEX {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val reference_genome
    shell:
    '''
    #!/usr/bin/env bash
    MEM=$(echo !{task.memory} | cut -d' ' -f1)
    REF=!{reference_genome}

    echo 'Index the reference genome for GATK4.'
    gatk --java-options "-Xmx${MEM}g" \
        CreateSequenceDictionary \
        -R ${REF} \

    echo "Output:"
    echo "  (1/1) {reference_genome_no_ext}.dict"
    '''
}

process JULIA_INSTALL_PACKAGES {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    echo 'using Pkg; Pkg.add(["StatsBase", "ProgressMeter", "DataFrames", "Plots"])' > install_julia_packages.jl
    julia install_julia_packages.jl
    rm install_julia_packages.jl
    '''
}

process INSTALL_POOLGEN {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
    output:
        val 0
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}
    cargo install --git https://github.com/jeffersonfparil/poolgen.git
    '''
}

workflow {
    GENOME_FIX(params.reference_genome) | \
    (GENOME_BWA_INDEX & GENOME_SAMTOOLS_INDEX & GENOME_GATK4_INDEX)
    JULIA_INSTALL_PACKAGES(params.dir_reads)
}