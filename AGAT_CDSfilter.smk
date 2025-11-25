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
    "N_forsteri",
    "test"
]


rule target:
    input:
        expand("results/tiberius/agat/{genome}.qc.yaml", genome=input_genomes),


rule summerize_qc:
    input:
        gtf="results/tiberius/agat/{genome}.qc.gff",
    output:
        yaml="results/tiberius/agat/{genome}.qc.yaml",
    resources:
        mem="32G",
        runtime=60,
    log:
        "logs/agat/{genome}.qc.log",
    container:
        agat
    shell:
        "agat_sp_statistics.pl "
        "--yaml "
        "--gff {input.gtf} "
        "--output results/tiberius/agat/{wildcards.genome}.qc "
        "&>> {log}"


rule report_CDS_phases:
    input:
        gff="results/tiberius/agat/{genome}.qc.gff.temp2",
        fasta="data/genomes/{genome}.fasta",
    output:
        gff="results/tiberius/agat/{genome}.qc.gff",
    resources:
        mem="16G",
        runtime=30,
    log:
        "logs/agat/{genome}.qc.log",
    container:
        agat
    shell:
        "agat_sp_fix_cds_phases.pl "
        "--gff {input.gff} "
        "--fasta {input.fasta} "
        "--output {output.gff} "
        "&>> {log}"


rule flag_premature_stop_codon:
    input:
        gff="results/tiberius/agat/{genome}.qc.gff.temp1",
        fasta="data/genomes/{genome}.fasta",
    output:
        gff="results/tiberius/agat/{genome}.qc.gff.temp2",
    resources:
        mem="16G",
        runtime=30,
    log:
        "logs/agat/{genome}.qc.log",
    container:
        agat
    shell:
        "agat_sp_flag_premature_stop_codons.pl "
        "--gff {input.gff} "
        "--fasta {input.fasta} "
        "--output {output.gff} "
        "&>> {log}"


rule filter_incomplete_CDS:
    input:
        gtf="results/tiberius/{genome}.gff",
        fasta="data/genomes/{genome}.fasta",
    output:
        gff="results/tiberius/agat/{genome}.qc.gff.temp1",
    resources:
        mem="16G",
        runtime=30,
    log:
        "logs/agat/{genome}.qc.log",
    container:
        agat
    shell:
        "agat_sp_filter_incomplete_gene_coding_models.pl "
        "--gff {input.gtf} "
        "--fasta {input.fasta} "
        "--add_flag "
        "--output {output.gff} "
        "&>> {log}"
