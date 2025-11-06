#!/usr/bin/env python3
def get_busco_db(wildcards):
    """Get the appropriate BUSCO database for each genome"""
    return f"data/busco_db/db/{genome_config[wildcards.genome]['busco_db']}"


def get_lineage_name(wildcards):
    """Get the lineage name for the JSON filename"""
    return genome_config[wildcards.genome]["busco_db"]


# containers
busco = "docker://ezlabgva/busco:v6.0.0_cv1"

# config
input_genomes = [
    "N_forsteri.8",
]

genome_config = {
    "A_magna": {"busco_db": "eukaryota_odb10"},
    "E_pictum": {"busco_db": "eukaryota_odb10"},
    "R_gram": {
        "busco_db": "helotiales_odb10",
    },
    "X_john": {"busco_db": "liliopsida_odb10"},
    "T_triandra": {"busco_db": "poales_odb10"},
    "H_bino": {"busco_db": "eukaryota_odb10"},
    "P_vit": {"busco_db": "eukaryota_odb10"},
    "P_halo": {"busco_db": "actinopterygii_odb10"},
    "N_erebi": {"busco_db": "actinopterygii_odb10"},
    "N_cryptoides": {"busco_db": "hymenoptera_odb10"},
    "N_forsteri.8": {"busco_db": "vertebrata_odb10"},
}

# db_path = "data/funannotate_db"


rule target:
    input:
        expand("results/tiberius/busco/{genome}/{genome}.json", genome=input_genomes),


rule busco:
    input:
        gtf="results/tiberius/{genome}.gtf",
    output:
        json="results/tiberius/busco/{genome}/{genome}.json",
    params:
        busco_db=get_busco_db,
        lineage=get_lineage_name,
        outdir="results/tiberius/busco/{genome}",
    resources:
        mem="32G",
        runtime=60,
    log:
        "logs/busco/{genome}.log",
    container:
        busco
    shell:
        "busco "
        "-i {input.gtf} "
        "-o {params.outdir} "
        "-l {params.busco_db} "
        "-m protein "
        "--force "
        "&> {log}; "
        "cp {params.outdir}/short_summary.specific.{params.lineage}.{wildcards.genome}.json {output.json}"
