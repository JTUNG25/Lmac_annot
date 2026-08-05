#!/usr/bin/env python3

tiberius = "/home/uqctung/containers/tiberius_latest.sif"


input_genomes = [
    "H_bac",
]


rule target:
    input:
        expand("results/tiberius/{genome}.gff3.gz", genome=input_genomes),
        expand("results/tiberius/{genome}.gtf.gz", genome=input_genomes),


rule compress_tiberius_output:
    input:
        gff3="results/tiberius/{genome}.gff3",
        gtf="results/tiberius/{genome}.gtf",
    output:
        gff3_gz="results/tiberius/{genome}.gff3.gz",
        gtf_gz="results/tiberius/{genome}.gtf.gz",
    log:
        "logs/tiberius/compressed_results/{genome}.log",
    resources:
        mem_mb=4000,
        runtime=20,
    shell:
        "gzip -k {input.gff3} {input.gtf} &> {log}"


rule tiberius:
    input:
        fasta="data/genomes/{genome}.fasta",
    output:
        gtf="results/tiberius/{genome}.gtf",
        gff3="results/tiberius/{genome}.gff3",
    log:
        "logs/tiberius/{genome}.log",
    container:
        tiberius
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
        "tiberius.py "
        "--genome {input.fasta} "
        "--model_cfg {params.model_cfg} "
        "--out {output.gtf} {output.gff3} "
        "--batch_size {params.batch_size} "
        "&> {log}"
