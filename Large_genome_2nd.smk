#!/usr/bin/env python3

# containers
tiberius = "docker://larsgabriel23/tiberius@sha256:796a9de5fdef73dd9148360dd22af0858c2fe8f0adc45ecaeda925ea4d4105d3"
bbmap = "docker://quay.io/biocontainers/bbmap:39.37--he5f24ec_0"  # new version for
# config
input_genomes = [
    "N_forsteri",
]

no_of_splits = 2
splits = [str(x) for x in range(0, no_of_splits, 1)]


wildcard_constraints:
    i="|".join(splits),
    genome="|".join(input_genomes),


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


rule gather_annotation:
    input:
        gtf=expand("results/tiberius/{{genome}}.20.{i}.gtf", i=splits),
    output:
        gtf="results/tiberius/{genome}.gtf",
    shell:
        "cat {input.gtf} > {output.gtf}"


rule tiberius:
    input:
        fasta="results/{genome}/partition/genome.20.{i}.fa",
        model="data/tiberius_weights_v2",
    output:
        gtf="results/tiberius/{genome}.20.{i}.gtf",
    params:
        #seq_len=259992,
        batch_size=8,
    resources:
        mem="360G",
        runtime=240,
        gpu=1,
        partitionFlag="--partition=gpu-a100",
        exclusive="--exclusive",
    log:
        "logs/tiberius/{genome}.20.{i}.log",
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


rule partition:
    input:
        "results/{genome}/partition/genome.20.fa",
    output:
        expand("results/{{genome}}/partition/genome.20.{i}.fa", i=splits),
    log:
        "logs/partition/{genome}.log",
    threads: 1
    params:
        pattern="results/{genome}/partition/genome.%.fa",
        ways=no_of_splits,
    resources:
        runtime=10,
        mem_mb=int(128e3),
    container:
        bbmap
    shell:
        "partition.sh -Xmx{resources.mem_mb}m "
        "bp=t "
        "ways={params.ways} "
        "in={input} "
        "out={params.pattern} 2>{log}"
