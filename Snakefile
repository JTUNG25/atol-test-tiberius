#!/usr/bin/env python3


rule target:
    input:
        "results/tiberius/RGram.gtf",


rule tiberius:
    input:
        fasta="data/{genome}.fasta",
        model="data/tiberius_weights_v2",
    output:
        gtf="results/tiberius/{genome}.gtf",
    params:
        seq_len=259992,
    threads: 1
    resources:
        mem="32G",
        runtime=10,
        gpu=1,
        partitionFlag="--partition=gpu-a100",
    log:
        "logs/tiberius/{genome}.log",
    container:
        "docker://quay.io/biocontainers/tiberius:1.1.6--pyhdfd78af_0"
    shell:
        # FIXME. python package doesn't get installed in biocontainer. Models
        # don't get shipped either. Provide the model weights (not config).
        # Find the weights URL in the config and download it manually.
        "python /usr/local/lib/python3.12/site-packages/tiberius/main.py "
        "--genome {input.fasta} "
        "--model {input.model} "
        "--out {output.gtf} "
        "--seq_len {params.seq_len} "
        "&> {log}"
