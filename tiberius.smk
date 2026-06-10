#!/usr/bin/env python3

# containers
tiberius = "docker://larsgabriel23/tiberius@sha256:c35ac0b456ee95df521e19abb062329fc8e39997723196172e10ae2c345f41e3"  # Nov 2025 updated container
bbmap = "docker://quay.io/biocontainers/bbmap:39.37--he5f24ec_0"  # new version for bp=t
# config
input_genomes = [
    "L_mac",
]


rule target:
    input:
        expand("results/tiberius/{genome}.gtf.gz", genome=input_genomes),


rule compress_tiberius_output:
    input:
        gtf="results/tiberius/{genome}.gtf",
    output:
        gtf_gz="results/tiberius/{genome}.gtf.gz",
    resources:
        mem="4G",
        runtime=20,
    log:
        "logs/tiberius/compressed_results/{genome}.log",
    container:
        tiberius
    shell:
        "gzip -k {input.gtf}"


rule tiberius:
    input:
        fasta="data/genomes/{genome}.fasta",
        model="data/tiberius_weights_v2",
    output:
        gtf="results/tiberius/{genome}.gtf",
    params:
        #seq_len=259992,
        batch_size=8,
    resources:
        mem="512G",
        runtime=1440,
        gpu=1,
        partitionFlag="--partition=gpu-h100",
        exclusive="--exclusive",
    benchmark:
        "benchmarks/tiberius/{genome}.tsv",
    log:
        "logs/tiberius/{genome}.log",
    container:
        # "docker://quay.io/biocontainers/tiberius:1.1.6--pyhdfd78af_0" FIXME.
        # The biocontainer tensorflow doesn't work, but the dev container
        # isn't versioned.
        tiberius
    shell:
        # FIXME. python package doesn't get installed in biocontainer. Models
        # don't get shipped either. Provide the model weights (not config).
        # "https://bioinf.uni-greifswald.de/bioinf/tiberius/models/tiberius_weights_v2.tar.gz"
        # Find the weights URL in the config and download it manually. This
        # needs to be checked for the dev container.
        "nvidia-smi && "
        "tiberius.py "
        "--genome {input.fasta} "
        "--model {input.model} "
        "--out {output.gtf} "
        "--batch_size {params.batch_size} "
        "&> {log}"
        #"--seq_len {params.seq_len} "

rule reformat:
    input:
        "data/genomes/{genome}.fasta",
    output:
        temp("results/{genome}/reformat/genome.fasta"),
    log:
        "logs/reformat/{genome}.log",
    threads: 1
    resources:
        runtime=10,
        mem_mb=int(32e3),
    container:
        bbmap
    shell:
        "reformat.sh -Xmx{resources.mem_mb}m "
        "in={input} out={output} "
        "addunderscore=t "
        "2>{log}"