#!/usr/bin/env snakemake -s
# ══════════════════════════════════════════════════════════════════════════════
# Tiberius Functional Annotation Pipeline (Dual Samples)
# AB INITIO vs EVIDENCE-BASED comparison
# eggNOG-mapper + DIAMOND BLAST → Combined Functional Report
# ══════════════════════════════════════════════════════════════════════════════

import os

# ── Container Images ──────────────────────────────────────────────────────────
eggnog = "/QRISdata/Q9140/lmac/annot_lmac/sifs/eggnog.sif"
diamond = "/QRISdata/Q9140/lmac/annot_lmac/sifs/diamond.sif"
biopython = "/QRISdata/Q9140/lmac/annot_lmac/sifs/biopython.sif"

# ── Sample Configuration ──────────────────────────────────────────────────────
SAMPLES = {
    "Lmac_D5": {
        "proteins": "/QRISdata/Q9140/lmac/annot_lmac/results/tiberius_evidence/Lmac_D5_proteins.fa",
        "gff": "/QRISdata/Q9140/lmac/annot_lmac/results/tiberius_evidence/Lmac_D5.gff3",
        "label": "Evidence-Based",
    },
    "Lmac_D5_ab": {
        "proteins": "/QRISdata/Q9140/lmac/annot_lmac/results/tiberius_evidence/Lmac_D5_ab_proteins.fa",
        "gff": "/QRISdata/Q9140/lmac/annot_lmac/results/tiberius_evidence/Lmac_D5_ab.gff3",
        "label": "AB Initio",
    },
}

# ── Databases ─────────────────────────────────────────────────────────────────
EGGNOG_DB = "/QRISdata/Q9140/lmac/annot_lmac/databases/eggnog"
DIAMOND_DB = "/QRISdata/Q9140/lmac/annot_lmac/databases/uniprot_sprot.dmnd"

EVALUE_CUTOFF = 1e-5
EGGNOG_TAX_SCOPE = 4751  # Fungi (NCBI taxon ID)


# ══════════════════════════════════════════════════════════════════════════════
# TARGET RULES
# ══════════════════════════════════════════════════════════════════════════════

rule target:
    input:
        expand("results/annotation/{sample}_functional_annotation.tsv", sample=SAMPLES.keys()),
        expand("results/annotation/{sample}_functional_report.xlsx", sample=SAMPLES.keys()),
        expand("results/stats/{sample}_annotation_summary.txt", sample=SAMPLES.keys()),


# ══════════════════════════════════════════════════════════════════════════════
# ANNOTATION RULES
# ══════════════════════════════════════════════════════════════════════════════

rule eggnog_annotation:
    input:
        proteins=lambda wildcards: SAMPLES[wildcards.sample]["proteins"],
    output:
        annotations="results/eggnog/{sample}.emapper.annotations",
        hits="results/eggnog/{sample}.emapper.hits",
        seed_orthologs="results/eggnog/{sample}.emapper.seed_orthologs",
    log:
        "logs/eggnog_{sample}.log",
    container:
        eggnog
    threads: 12
    resources:
        mem_mb=32000,
        runtime=240,
    shell:
        """
        mkdir -p results/eggnog
        
        export TMPDIR=/QRISdata/Q9140/lmac/annot_lmac/results/tiberius_evidence/tmp
        mkdir -p $TMPDIR
        
        emapper.py \
            -i {input.proteins} \
            --data_dir {EGGNOG_DB} \
            --output /QRISdata/Q9140/lmac/annot_lmac/results/tiberius_evidence/results/eggnog/{wildcards.sample} \
            -m diamond \
            --tax_scope {EGGNOG_TAX_SCOPE} \
            --go_evidence non-electronic \
            --target_orthologs all \
            --seed_ortholog_evalue 0.001 \
            --seed_ortholog_score 60 \
            --override \
            --cpu {threads} \
            >{log} 2>&1
        """


rule diamond_blast:
    input:
        proteins=lambda wildcards: SAMPLES[wildcards.sample]["proteins"],
        database=DIAMOND_DB,
    output:
        results="results/diamond/{sample}_blast.tsv",
    log:
        "logs/diamond_{sample}.log",
    container:
        diamond
    threads: 12
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        """
        mkdir -p results/diamond
        
        diamond blastp \
            -d {input.database} \
            -q {input.proteins} \
            -o {output.results} \
            -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
            --max-target-seqs 5 \
            --evalue {EVALUE_CUTOFF} \
            --threads {threads} \
            >{log} 2>&1
        """


rule combine_annotations:
    input:
        eggnog="results/eggnog/{sample}.emapper.annotations",
        diamond="results/diamond/{sample}_blast.tsv",
        gff=lambda wildcards: SAMPLES[wildcards.sample]["gff"],
    output:
        combined="results/annotation/{sample}_functional_annotation.tsv",
    log:
        "logs/combine_annotations_{sample}.log",
    container:
        biopython
    threads: 4
    resources:
        mem_mb=8000,
        runtime=30,
    script:
        "scripts/combine_annotations_tiberius.py"


rule generate_report:
    input:
        combined="results/annotation/{sample}_functional_annotation.tsv",
    output:
        report="results/annotation/{sample}_functional_report.xlsx",
        summary="results/stats/{sample}_annotation_summary.txt",
    log:
        "logs/generate_report_{sample}.log",
    container:
        biopython
    threads: 2
    resources:
        mem_mb=4000,
        runtime=30,
    script:
        "scripts/generate_report_tiberius.py"
