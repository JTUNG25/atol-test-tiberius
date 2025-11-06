#!/usr/bin/env python3

# containers
agat = "docker://quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0"

# config
input_genomes = [
    "A_magna",
    "E_pictum",
    "R_gram",
    "X_john",
    "T_triandra",
    "H_bino",
    "P_vit",
    "P_halo",
    "N_erebi",
    "N_cryptoides",
    "N_forsteri.8",
]


rule target:
    input:
        expand("results/tiberius/agat/{genome}.yaml", genome=input_genomes),

rule agat:
    input:
        gtf="results/tiberius/{genome}.gtf",
    output:
        txt="results/tiberius/agat/{genome}.yaml",
    resources:
        mem="32G",
        runtime=60,
    log:
        "logs/agat/{genome}.log",
    container:
        agat
    shell:
        "agat_sp_statistics.pl "
        "--yaml "
        "--gff {input.gtf} "
        "--output {output.txt} "
        "&> {log}"