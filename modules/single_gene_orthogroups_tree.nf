//////////////////////////////////
// SINGLE-GENE ORTHOGROUPS TREE //
//////////////////////////////////
// Assumed extension names:
//  - genome: *.fna
//  - annotation: *.gff
//  - coding DNA: *.cds
//  - proteome: *.faa

process IDENTIFY_SINGLE_GENE_ORTHOGROUPS {
    label "HIGH_MEM_HIGH_CPU"
    input:
        val dir
        val dates
    output:
        val dir
        val dates
    shell:
    '''
    #!/usr/bin/env bash
    cd !{dir}

    
    trimmomatic \
        PE \
        -threads 32 \
        reads/test_R1.fastq.gz        reads/test_R2.fastq.gz \
        reads/test_R1_paired.fastq.gz reads/test_R1_unpaired.fastq.gz \
        reads/test_R2_paired.fastq.gz reads/test_R2_unpaired.fastq.gz \
        ILLUMINACLIP:IDT_for_Illumina_TruSeq_UD_and_CD_indexes.fa:2:30:10:2:True \
        LEADING:3 \
        TRAILING:3 \
        MINLEN:36

    echo "Cleanup"
    rm single_gene_list.grep 
    rm parallel_extract_single_gene_orthogroups.sh
    
    echo "Output:"
    echo "  (1/2) {ORTHONAME}.fasta"
    echo "  (2/2) {ORTHONAME}-{SPECIES}.fasta"
    '''
}

workflow {
    IDENTIFY_SINGLE_GENE_ORTHOGROUPS(params.dir, params.dates) | \
        ALIGN_SINGLE_GENE_ORTHOGROUPS | \
        BUILD_TREE
}
