#!/usr/bin/env python3

#"docker://larsgabriel23/tiberius@sha256:c35ac0b456ee95df521e19abb062329fc8e39997723196172e10ae2c345f41e3" older version that need update


input_genomes = [
    "Lmac_D5",
]


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
        fasta="data/genomes/{genome}.fasta",
    output:
        gtf="results/tiberius/{genome}.gtf",
    log:
        "logs/tiberius/{genome}.log",
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
        "source tiberius_venv/bin/activate && "
        "python Tiberius/tiberius.py "
        "--genome {input.fasta} "
        "--model_cfg Tiberius/model_cfg/{params.model_cfg}.yaml "
        "--out {output.gtf} "
        "--batch_size {params.batch_size} "
        "&> {log}"