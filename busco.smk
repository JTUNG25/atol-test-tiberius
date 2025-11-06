#!/usr/bin/env python3

# containers
agat = "docker:ezlabgva/busco:v6.0.0_cv1"

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

#db_path = "data/funannotate_db"


rule target:
    input:
        expand("results/tiberius/busco/{genome}.json", genome=input_genomes),


rule agat:
    input:
        gtf="results/tiberius/{genome}.gtf",
    output:
        json="results/tiberius/busco/{genome}.json",
    params:
        busco_db="data/busco_db/db/vertebrata_odb10"
    resources:
        mem="32G",
        runtime=60,
    log:
        "logs/busco/{genome}.log",
    container:
        agat
    shell:
        "--i {input.gtf} "
        "--o {output.json} "
        "--lineage_dataset {params.busco_db} "
        "--mode protein "
        "&> {log}"
