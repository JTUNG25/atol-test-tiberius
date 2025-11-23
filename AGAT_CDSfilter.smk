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
]


rule target:
    input:
        expand("results/tiberius/agat/{genome}.qc.tsv", genome=input_genomes),

rule summerize_qc:
    input:
        gff="results/tiberius/agat/{genome}.qc.gff",
    output:
        tsv="results/tiberius/agat/{genome}.qc.tsv",
    resources:
        mem="16G",
        runtime=60,
    log:
        "logs/agat/{genome}.qc.log",
    shell:
        "echo -e 'genome\\ttotal_mRNA\\tmissing_start\\tpct_start\\tmissing_stop\\tpct_stop\\tboth_missing\\tinternal_stop\\tpct_internal\\tcds_not_div_by_3\\tpct_not_div_by_3' > 2>&1 | tee {log} && "
        "total=$(awk '$3 ~ /mRNA|transcript/' {input.gff} | wc -l); [[ $total -eq 0 ]] && total=1; "
        "missing_start=$(grep -c 'incomplete=1' {input.gff} || echo 0); "
        "missing_stop=$(grep -c 'incomplete=2' {input.gff} || echo 0); "
        "both_missing=$(grep -c 'incomplete=3' {input.gff} || echo 0); "
        "internal_stop=$(grep -c 'pseudo=1' {input.gff} || echo 0); "
        "not_div_by_3=$(awk '$3==\"CDS\" {{len=$5-$4+1; if(len%3!=0) print $9}}' {input.gff} | grep -o 'Parent=[^;]*' | sort -u | wc -l); "
        "pct_start=$(awk \"BEGIN {{printf \\\"%.2f\\\", 100*$missing_start/$total}}\"); "
        "pct_stop=$(awk \"BEGIN {{printf \\\"%.2f\\\", 100*$missing_stop/$total}}\"); "
        "pct_both=$(awk \"BEGIN {{printf \\\"%.2f\\\", 100*$both_missing/$total}}\"); "
        "pct_int=$(awk \"BEGIN {{printf \\\"%.2f\\\", 100*$internal_stop/$total}}\"); "
        "pct_div3=$(awk \"BEGIN {{printf \\\"%.2f\\\", 100*$not_div_by_3/$total}}\"); "
        "echo -e \"{wildcards.genome}\\t$total\\t$missing_start\\t$pct_start\\t$missing_stop\\t$pct_stop\\t$both_missing\\t$internal_stop\\t$pct_int\\t$not_div_by_3\\t$pct_div3\" >> {output.tsv} 2>&1 | tee -a {log}"


rule agat_filter_incomplete_CDS:
    input:
        gtf="results/tiberius/{genome}.gtf",
        fasta="data/genomes/{genome}.fasta",
    output:
        gff="results/tiberius/agat/{genome}.qc.gff",
    resources:
        mem="32G",
        runtime=120,
    log:
        "logs/agat/{genome}.log",
    container:
        agat
    shell:
        "agat_sp_filter_incomplete_gene_coding_models.pl "
        "--gff {input.gtf} "
        "--fasta {input.fasta} "
        "--add_flag "
        "--output {output.gff}.temp1 "
        "&>> {log} && "
        "agat_sp_flag_premature_stop_codons.pl "
        "--gff {output.gff}.temp1 "
        "--fasta {input.fasta} "
        "--output {output.gff}.temp2 "
        "&>> {log} && "
        # 3. fix/report CDS phases (detects length % 3 != 0)
        "agat_sp_fix_cds_phases.pl "
        "--gff {output.gff}.temp2 "
        "--fasta {input.fasta} "
        "--output {output.gff} "
        "&>> {log} && "
        "rm -f {output.gff}.temp1 {output.gff}.temp2"