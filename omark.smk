#!/usr/bin/env python3

# containers
omamer = "docker://quay.io/biocontainers/omamer:2.1.2--pyhdfd78af_0"
omark = "/scratch/user/uqctung/apptainer_cache/d50221fcdbdecb24f5454545f68e896a.simg"
#"docker://quay.io/biocontainers/omark:0.4.1--pyh7e72e81_0"

GENOME = "Lmac_D5_CLEANED"  
PROTEINS = "/QRISdata/Q9140/lmac/annot_lmac/results/tiberius_evidence/tiberius_train_proteins.fa"
OMA_DB = "data/omark/Saccharomyceta.h5 "


rule target:
    input:
        f"results/omark/{GENOME}/{GENOME}.sum",

rule omamer_search:
    input:
        proteins=PROTEINS,
        db=OMA_DB,
    output:
        omamer_out=f"results/omark/{GENOME}/{GENOME}.omamer",
    log:
        f"logs/omark/{GENOME}.omamer.log",
    container:
        omamer
    threads: 8
    resources:
        mem_mb=64000,
        runtime=120,
    shell:
        "omamer search "
        "--db {input.db} "
        "--query {input.proteins} "
        "--out {output.omamer_out} "
        "--nthreads {threads} "
        "&> {log}"


rule omark:
    input:
        omamer_out=f"results/omark/{GENOME}/{GENOME}.omamer",
        db=OMA_DB,
        proteins=PROTEINS,
    output:
        f"results/omark/{GENOME}/{GENOME}.sum",
    log:
        f"logs/omark/{GENOME}.omark.log",
    container:
        omark
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        outdir=f"results/omark/{GENOME}",
        ete_db="data/omark/ete3/taxa.sqlite",
    shell:
        "omark "
        "-f {input.omamer_out} "
        "-d {input.db} "
        "-of {input.proteins} "
        "-o {params.outdir} "
        "-e {params.ete_db} "
        "&> {log}"