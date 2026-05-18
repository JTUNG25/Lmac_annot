#!/usr/bin/env python3

eggnog = "/QRISdata/Q9140/lmac/annot_lmac/sifs/eggnog.sif"
interproscan = "/QRISdata/Q9140/lmac/annot_lmac/sifs/interproscan.sif"
diamond = "/QRISdata/Q9140/lmac/annot_lmac/sifs/diamond.sif"
biopython = "/QRISdata/Q9140/lmac/annot_lmac/sifs/biopython.sif"

GENOME = "data/genome/JN3.fasta"
ANNOTATION = "data/genome/JN3.gff"
MUTANTS = {
    # ── AGO mutants ──────────────────────────────────────────────────────────
    "ago1_silencing": "data/target_lists/ago1_silencing.txt",
    "ago1_derepression": "data/target_lists/ago1_derepression.txt",
    "ago13_silencing": "data/target_lists/ago13_silencing.txt",
    "ago13_derepression": "data/target_lists/ago13_derepression.txt",
    "ago3_silencing": "data/target_lists/ago3_silencing.txt",
    "ago3_derepression": "data/target_lists/ago3_derepression.txt",
    # ── DCL mutants ──────────────────────────────────────────────────────────
    "dcl1_silencing": "data/target_lists/dcl1_silencing.txt",
    "dcl1_derepression": "data/target_lists/dcl1_derepression.txt",
    "dcl2_derepression": "data/target_lists/dcl2_derepression.txt",
    # dcl2_silencing: not generated — no Silencing pairs in strict consensus
    # ── RDRP mutants ─────────────────────────────────────────────────────────
    "rdrp1_silencing": "data/target_lists/rdrp1_silencing.txt",
    "rdrp1_derepression": "data/target_lists/rdrp1_derepression.txt",
    "rdrp2_silencing": "data/target_lists/rdrp2_silencing.txt",
    "rdrp2_derepression": "data/target_lists/rdrp2_derepression.txt",
    # rdrp12: not generated — all strict consensus pairs were Basal Silencing
}

# Pipeline settings
RUN_INTERPROSCAN = (
    False  # Set to True if you want domain analysis (adds 3+ hours per mutant)
)
MAX_EVALUE = 1e-5

# Get all mutant names
mutant_names = list(MUTANTS.keys())


rule target:
    input:
        expand("results/annotation/{mutant}/final_annotation.xlsx", mutant=mutant_names),
        expand("results/annotation/{mutant}/readable_summary.csv", mutant=mutant_names),
        "results/summary/all_mutants_annotation_summary.xlsx",


rule extract_sequences:
    input:
        genome=GENOME,
        annotation=ANNOTATION,
        gene_list=lambda wildcards: MUTANTS[wildcards.mutant],
    output:
        proteins="results/sequences/{mutant}/proteins.fasta",
        cds="results/sequences/{mutant}/cds.fasta",
        summary="results/sequences/{mutant}/extraction_summary.csv",
    log:
        "logs/extract_sequences/{mutant}.log",
    container:
        biopython
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
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
        proteins="results/sequences/{mutant}/proteins.fasta",
        data_dir="databases/eggnog",
    output:
        annotations="results/eggnog/{mutant}/eggnog.emapper.annotations",
        hits="results/eggnog/{mutant}/eggnog.emapper.hits",
        seed_orthologs="results/eggnog/{mutant}/eggnog.emapper.seed_orthologs",
    log:
        "logs/eggnog/{mutant}.log",
    container:
        eggnog
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        """
        mkdir -p results/eggnog/{wildcards.mutant}

        emapper.py \
        -i {input.proteins} \
        --data_dir /QRISdata/Q9140/lmac/annot_lmac/databases/eggnog \
        --output /QRISdata/Q9140/lmac/annot_lmac/results/eggnog/{wildcards.mutant}/eggnog \
        -m diamond \
        --tax_scope 4751 \
        --go_evidence non-electronic \
        --target_orthologs all \
        --seed_ortholog_evalue 0.001 \
        --seed_ortholog_score 60 \
        --override \
        --cpu {threads} \
        > {log} 2>&1
        """


rule interproscan:
    input:
        proteins="results/sequences/{mutant}/proteins.fasta",
    output:
        results=(
            "results/interproscan/{mutant}/interproscan_results.tsv"
            if RUN_INTERPROSCAN
            else []
        ),
    log:
        "logs/interproscan/{mutant}.log",
    container:
        interproscan
    threads: 8
    resources:
        mem_mb=32000,
        runtime=360,
    shell:
        """
        mkdir -p results/interproscan/{wildcards.mutant}

        interproscan.sh \
            -i {input.proteins} \
            -f tsv \
            -o {output.results} \
            --goterms \
            --pathways \
            --cpu {threads} \
            > {log} 2>&1
        """
        if RUN_INTERPROSCAN
        else "echo 'InterProScan skipped' > {log}"


rule diamond_blast:
    input:
        proteins="results/sequences/{mutant}/proteins.fasta",
        database="databases/uniprot_sprot.dmnd",
    output:
        results="results/diamond/{mutant}/blast_results.tsv",
    log:
        "logs/diamond/{mutant}.log",
    container:
        diamond
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
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
        interproscan=(
            "results/interproscan/{mutant}/interproscan_results.tsv"
            if RUN_INTERPROSCAN
            else []
        ),
        gene_list=lambda wildcards: MUTANTS[wildcards.mutant],
    output:
        combined="results/annotation/{mutant}/combined_annotations.xlsx",
        summary="results/annotation/{mutant}/annotation_summary.csv",
    log:
        "logs/combine_annotations/{mutant}.log",
    container:
        biopython
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    script:
        "scripts/combine_annotations.py"


rule final_annotation:
    input:
        combined="results/annotation/{mutant}/combined_annotations.xlsx",
        gene_list=lambda wildcards: MUTANTS[wildcards.mutant],
    output:
        final="results/annotation/{mutant}/final_annotation.xlsx",
        readable="results/annotation/{mutant}/readable_summary.csv",
    log:
        "logs/final_annotation/{mutant}.log",
    container:
        biopython
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15,
    script:
        "scripts/create_final_report.py"


rule summarize_all_mutants:
    input:
        annotations=expand(
            "results/annotation/{mutant}/final_annotation.xlsx", mutant=mutant_names
        ),
    output:
        summary="results/summary/all_mutants_annotation_summary.xlsx",
        comparison="results/summary/mutant_comparison.csv",
    log:
        "logs/summarize_all_mutants.log",
    container:
        biopython
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    run:
        import pandas as pd
        from pathlib import Path

        Path("results/summary").mkdir(exist_ok=True)

        all_data = []
        mutant_summaries = []

        for i, mutant in enumerate(mutant_names):
            try:
                df = pd.read_excel(input.annotations[i])
                df["Mutant"] = mutant
                df["Mutant_Type"] = "_".join(
                    mutant.split("_")[:-1]
                )  # ago1, dcl2, etc.
                df["Effect_Type"] = mutant.split("_")[
                    -1
                ]  # silencing / derepression

                all_data.append(df)

                summary_row = {
                    "Mutant": mutant,
                    "Mutant_Type": "_".join(mutant.split("_")[:-1]),
                    "Effect_Type": mutant.split("_")[-1],
                    "Total_Genes": len(df),
                    "With_Function": (
                        df["Predicted_Function"].notna().sum()
                        if "Predicted_Function" in df.columns
                        else 0
                    ),
                    "With_GO_Terms": (
                        df["GO_Terms"].notna().sum()
                        if "GO_Terms" in df.columns
                        else 0
                    ),
                    "High_Confidence": (
                        (df["Confidence"] == "High").sum()
                        if "Confidence" in df.columns
                        else 0
                    ),
                }
                mutant_summaries.append(summary_row)

            except Exception as e:
                print(f"Error processing {mutant}: {e}")

        if all_data:
            combined_df = pd.concat(all_data, ignore_index=True)
            summary_df = pd.DataFrame(mutant_summaries)

            with pd.ExcelWriter(output.summary, engine="openpyxl") as writer:
                combined_df.to_excel(
                    writer, sheet_name="All_Annotations", index=False
                )
                summary_df.to_excel(
                    writer, sheet_name="Mutant_Summary", index=False
                )

                if len(summary_df) > 0:
                    pivot_genes = summary_df.pivot(
                        index="Mutant_Type",
                        columns="Effect_Type",
                        values="Total_Genes",
                    ).fillna(0)
                    pivot_annotated = summary_df.pivot(
                        index="Mutant_Type",
                        columns="Effect_Type",
                        values="With_Function",
                    ).fillna(0)

                    pivot_genes.to_excel(writer, sheet_name="Genes_by_Type")
                    pivot_annotated.to_excel(writer, sheet_name="Annotated_by_Type")

            summary_df.to_csv(output.comparison, index=False)

            print(
                f"Summary saved with {len(combined_df)} total annotations "
                f"across {len(mutant_summaries)} mutant conditions"
            )


rule setup_databases:
    output:
        uniprot_db="databases/uniprot_sprot.dmnd",
    log:
        "logs/setup_databases.log",
    container:
        diamond
    threads: 4
    resources:
        mem_mb=16000,
        runtime=240,
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
