#!/usr/bin/env python3

import glob
import os

# containers
tiberius = "docker://larsgabriel23/tiberius@sha256:c35ac0b456ee95df521e19abb062329fc8e39997723196172e10ae2c345f41e3"  # Nov 2025 updated container
# "docker://larsgabriel23/tiberius@sha256:796a9de5fdef73dd9148360dd22af0858c2fe8f0adc45ecaeda925ea4d4105d3" # Can't handle contig > 550Mbp
# "docker://larsgabriel23/tiberius:1.1.7"  # newer container

bbmap = "docker://quay.io/biocontainers/bbmap:39.37--he5f24ec_0"  # new version for bp=t
# config
input_contigs = [
    "chr1_1",
]


wildcard_constraints:
    contig="|".join(input_contigs),


# main


rule target:
    input:
        expand("results/tiberius/{contig}.gtf.gz", contig=input_contigs),


rule compress_tiberius_output:
    input:
        gtf="results/tiberius/{contig}.gtf",
    output:
        gtf_gz="results/tiberius/{contig}.gtf.gz",
    resources:
        mem="4G",
        runtime=20,
    log:
        "logs/tiberius/compressed_results/{contig}.log",
    container:
        tiberius
    shell:
        "gzip -k {input.gtf}"


rule tiberius:
    input:
        fasta="results/{contig}/reformat/contig.fna",
        model="data/tiberius_weights_v2",
    output:
        gtf="results/tiberius/{contig}.gtf",
    params:
        #seq_len=259992,
        batch_size=16,
    resources:
        mem="512G",
        runtime=180,
        gpu=1,
        partitionFlag="--partition=gpu-h100",
        exclusive="--exclusive",
    benchmark:
        "benchmarks/tiberius/{contig}.tsv",
    log:
        "logs/tiberius/{contig}.log",
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
        "--keep-incomplete "
        "&> {log}"


rule reformat:
    input:
        "data/genomes/{contig}.fna",
    output:
        temp("results/{contig}/reformat/contig.fna"),
    log:
        "logs/reformat/{contig}.log",
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
