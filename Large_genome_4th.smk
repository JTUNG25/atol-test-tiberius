#!/usr/bin/env python3

import glob
import os

# containers

tiberius = "docker://larsgabriel23/tiberius@sha256:796a9de5fdef73dd9148360dd22af0858c2fe8f0adc45ecaeda925ea4d4105d3"
bbmap = "docker://quay.io/biocontainers/bbmap:39.37--he5f24ec_0"  # new version for bp=t


# config

input_genomes = [
    "N_forsteri",
]

input_sequences = [str(i) for i in [22]]


wildcard_constraints:
    genome="|".join(input_genomes),
    sequence="|".join(input_sequences),


#############
# FUNCTIONS #
#############


def get_tiberius_output(wildcards):
    """Get all partition files for ALL sequences after checkpoint completes"""
    all_gtf_gzs = []
    for seq in input_sequences:
        checkpoint_output = checkpoints.partition_sequences.get(
            genome=wildcards.genome, sequence=seq
        ).output[0]

        file_pattern = os.path.join(checkpoint_output, f"genome.{seq}.shred.{{chunk}}.fa")
        chunks = glob_wildcards(file_pattern).chunk

        # error handling
        if not chunks:
            print(f"Warning: No chunks found for {seq} in {checkpoint_output}")
            print(f"Looking for pattern: {file_pattern}")
            # List what's actually there
            if os.path.exists(checkpoint_output):
                print(f"Files in directory: {os.listdir(checkpoint_output)}")

        gtf_gzs = expand(
            "results/tiberius/{genome}.genome.{sequence}.shred.{chunk}.gtf.gz",
            chunk=chunks,
            sequence=seq,
            genome=wildcards.genome,
        )
        all_gtf_gzs.extend(gtf_gzs)

    return all_gtf_gzs


##############
# MAIN RULES #
##############


rule target:
    input:
        expand(
            "results/{genome}/all_done.txt",
            genome=input_genomes,
        ),


rule collect_results:
    input:
        all_gtf_gzs=get_tiberius_output,
    output:
        "results/{genome}/all_done.txt",
    resources:
        mem="8G",
        runtime=10,
    shell:
        "touch {output}"


rule compress_tiberius_output:
    input:
        gtf="results/tiberius/{genome}.genome.{sequence}.shred.{chunk}.gtf",
    output:
        gtf_gz="results/tiberius/{genome}.genome.{sequence}.shred.{chunk}.gtf.gz",
    resources:
        mem="4G",
        runtime=10,
    log:
        "logs/tiberius/compressed_results/{genome}.{sequence}.{chunk}.log",
    shell:
        "gzip -c {input.gtf} > {output.gtf_gz} 2> {log}"


rule tiberius:
    input:
        fasta="results/{genome}/partition/partition2/{sequence}/genome.{sequence}.shred.{chunk}.fa",
        model="data/tiberius_weights_v2",
    output:
        gtf="results/tiberius/{genome}.genome.{sequence}.shred.{chunk}.gtf",
    params:
        #seq_len=259992,
        batch_size=8,
    resources:
        mem="256G",
        runtime=180,
        gpu=1,
        partitionFlag="--partition=gpu-a100",
        exclusive="--exclusive",
    log:
        "logs/tiberius/{genome}.{sequence}.{chunk}.log",
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


checkpoint partition_sequences:
    input:
        "results/{genome}/partition/genome.{sequence}.shred.fa",
    output:
        directory("results/{genome}/partition/partition2/{sequence}/"),
    params:
        pattern="results/{genome}/partition/partition2/{sequence}/genome.{sequence}.shred.%.fa",
    log:
        "logs/partition/{genome}.{sequence}.partition2.log",
    threads: 1
    resources:
        runtime=10,
        mem_mb=int(32e3),
    container:
        bbmap
    shell:
        "mkdir -p {output} && "
        'No_Seqs=$(grep -c "^>" {input}) && '
        "partition.sh -Xmx{resources.mem_mb}m "
        "in={input} "
        "out={params.pattern} "
        "ways=$No_Seqs "
        "2>{log}"


rule shred:
    input:
        "results/{genome}/partition/genome.{sequence}.fa",
    output:
        "results/{genome}/partition/genome.{sequence}.shred.fa",
    log:
        "logs/",
    threads: 1
    resources:
        runtime=10,
        mem_mb=int(32e3),
    container:
        bbmap
    shell:
        "shred.sh -Xmx{resources.mem_mb}m "
        "length=550000000 "
        "overlap=1000000 "
        "equal=f "
        "in={input} "
        "out={output} 2>{log}"
