#!/usr/bin/env python3

#"docker://larsgabriel23/tiberius@sha256:c35ac0b456ee95df521e19abb062329fc8e39997723196172e10ae2c345f41e3" older version that need update


input_genomes = [
    "Lmac_D5",
]

PROJECT_DIR = "/QRISdata/Q9140/lmac/annot_lmac"


rule target:
    input:
        expand("results/tiberius/{genome}.gtf.gz", genome=input_genomes),


rule compress_tiberius_output:
    input:
        gtf="results/tiberius/{genome}.gtf",
    output:
        gtf_gz="results/tiberius/{genome}.gtf.gz",
    log:
        "logs/tiberius/compressed_results/{genome}.log",
    resources:
        mem_mb=4000,
        runtime=20,
    shell:
        "gzip -k {input.gtf}"


rule tiberius:
    input:
        fasta=f"{PROJECT_DIR}/data/genomes/{{genome}}.fasta",
    output:
        gtf=f"{PROJECT_DIR}/results/tiberius/{{genome}}.gtf",
    log:
        f"{PROJECT_DIR}/logs/tiberius/{{genome}}.log",
    resources:
        mem_mb=500000,
        runtime=1440,
        partition="gpu_cuda",
        qos="gpu",
        gres="--gres=gpu:h100:1",
    params:
        batch_size=8,
        model_cfg="fungi",
    shell:
        "source {PROJECT_DIR}/tiberius_venv/bin/activate && "
        "python {PROJECT_DIR}/Tiberius/tiberius.py "
        "--genome {input.fasta} "
        "--model_cfg {PROJECT_DIR}/Tiberius/model_cfg/{params.model_cfg}.yaml "
        "--out {output.gtf} "
        "--batch_size {params.batch_size} "
        "&> {log}"