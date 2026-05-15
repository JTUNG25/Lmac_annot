#!/usr/bin/env python3

# Container definitions
eggnog = "docker://quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
interproscan = "docker://quay.io/biocontainers/interproscan:5.62_94.0--hec16e2b_1"
diamond = "docker://quay.io/biocontainers/diamond:2.1.8--h43eeafb_0"
biopython = "docker://quay.io/biocontainers/biopython:1.81"

# Input files - EDIT THESE PATHS
GENOME = "data/genome/JN3.fasta"
ANNOTATION = "data/genome/JN3_annotation.gff3"

# Mutant gene lists - EDIT THESE
MUTANTS = {
    "rnai_targets": "data/gene_lists/rnai_genes.txt",
}

# Pipeline settings
RUN_INTERPROSCAN = False  # Set to True for full analysis (adds 3+ hours)
MAX_EVALUE = 1e-5

# Get all mutant names
mutant_names = list(MUTANTS.keys())

rule target:
    input:
        expand("results/annotation/{mutant}/final_annotation.xlsx", mutant=mutant_names),
        expand("results/annotation/{mutant}/readable_summary.csv", mutant=mutant_names)

rule extract_sequences:
    input:
        genome=GENOME,
        annotation=ANNOTATION,
        gene_list=lambda wildcards: MUTANTS[wildcards.mutant]
    output:
        proteins="results/sequences/{mutant}/proteins.fasta",
        cds="results/sequences/{mutant}/cds.fasta",
        summary="results/sequences/{mutant}/extraction_summary.csv"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    container:
        biopython
    log:
        "logs/extract_sequences/{mutant}.log"
    shell:
        """
        mkdir -p results/sequences/{wildcards.mutant}
        
        python scripts/extract_sequences.py {input.genome} {input.annotation} \
            --genes {input.gene_list} \
            --output results/sequences/{wildcards.mutant} \
            > {log} 2>&1
        """

rule eggnog_annotation:
    input:
        proteins="results/sequences/{mutant}/proteins.fasta"
    output:
        annotations="results/eggnog/{mutant}/eggnog.emapper.annotations",
        hits="results/eggnog/{mutant}/eggnog.emapper.hits",
        seed_orthologs="results/eggnog/{mutant}/eggnog.emapper.seed_orthologs"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120,
    container:
        eggnog
    log:
        "logs/eggnog/{mutant}.log"
    shell:
        """
        mkdir -p results/eggnog/{wildcards.mutant}
        
        emapper.py \
            -i {input.proteins} \
            --output results/eggnog/{wildcards.mutant}/eggnog \
            --output_dir results/eggnog/{wildcards.mutant} \
            -m diamond \
            --tax_scope fungi \
            --go_evidence non-electronic \
            --target_orthologs all \
            --seed_ortholog_evalue 0.001 \
            --seed_ortholog_score 60 \
            --override \
            --cpu {threads} \
            > {log} 2>&1
        """

rule diamond_blast:
    input:
        proteins="results/sequences/{mutant}/proteins.fasta",
        database="databases/uniprot_sprot.dmnd"
    output:
        results="results/diamond/{mutant}/blast_results.tsv"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
    container:
        diamond
    log:
        "logs/diamond/{mutant}.log"
    shell:
        """
        mkdir -p results/diamond/{wildcards.mutant}
        
        diamond blastp \
            -d {input.database} \
            -q {input.proteins} \
            -o {output.results} \
            -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
            --max-target-seqs 5 \
            --evalue {MAX_EVALUE} \
            --threads {threads} \
            > {log} 2>&1
        """

rule combine_annotations:
    input:
        eggnog="results/eggnog/{mutant}/eggnog.emapper.annotations",
        diamond="results/diamond/{mutant}/blast_results.tsv",
        gene_list=lambda wildcards: MUTANTS[wildcards.mutant]
    output:
        combined="results/annotation/{mutant}/combined_annotations.xlsx",
        summary="results/annotation/{mutant}/annotation_summary.csv"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    container:
        biopython
    log:
        "logs/combine_annotations/{mutant}.log"
    script:
        "scripts/combine_annotations.py"

rule final_annotation:
    input:
        combined="results/annotation/{mutant}/combined_annotations.xlsx",
        gene_list=lambda wildcards: MUTANTS[wildcards.mutant]
    output:
        final="results/annotation/{mutant}/final_annotation.xlsx",
        readable="results/annotation/{mutant}/readable_summary.csv"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15,
    container:
        biopython
    log:
        "logs/final_annotation/{mutant}.log"
    script:
        "scripts/create_final_report.py"

rule setup_databases:
    output:
        uniprot_db="databases/uniprot_sprot.dmnd"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=240,
    container:
        diamond
    log:
        "logs/setup_databases.log"
    shell:
        """
        mkdir -p databases
        cd databases
        
        wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        gunzip uniprot_sprot.fasta.gz
        
        diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot
        rm uniprot_sprot.fasta
        
        echo "Database setup complete" > ../logs/setup_databases.log
        """
