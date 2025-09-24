#!/usr/bin/env python3
input_genomes = [
    "A_magna",
    "E_pictum",
    "R_gram",
    "X_john",
    "T_triandra",
    "H_bino",
    "P_vit",
    "P_halo",
    "N_erebi"
]


rule target:
    input:
        expand("results/tiberius/{genome}.gtf", genome=input_genomes),


rule tiberius:
    input:
        fasta="data/genome/{genome}.fasta",
        model="data/tiberius_weights_v2",
    output:
        gtf="results/tiberius/{genome}.gtf",
    params:
        seq_len=259992,
    resources:
        mem="256G",
        runtime=240,
        gpu=1,
        partitionFlag="--partition=gpu-a100",
    log:
        "logs/tiberius/{genome}.log",
    container:
        # "docker://quay.io/biocontainers/tiberius:1.1.6--pyhdfd78af_0" FIXME.
        # The biocontainer tensorflow doesn't work, but the dev container
        # isn't versioned.
        "docker://larsgabriel23/tiberius@sha256:796a9de5fdef73dd9148360dd22af0858c2fe8f0adc45ecaeda925ea4d4105d3"
    shell:
        # FIXME. python package doesn't get installed in biocontainer. Models
        # don't get shipped either. Provide the model weights (not config).
        # "https://bioinf.uni-greifswald.de/bioinf/tiberius/models/tiberius_weights_v2.tar.gz"
        # Find the weights URL in the config and download it manually. This
        # needs to be checked for the dev container.
        "tiberius.py "
        "--genome {input.fasta} "
        "--model {input.model} "
        "--out {output.gtf} "
        "--seq_len {params.seq_len} "
        "&> {log}"
