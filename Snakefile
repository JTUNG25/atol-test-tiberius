#!/usr/bin/env python3


rule target:
    input:
        "results/tiberius/RGram.gtf",


rule tiberius:
    input:
        fasta="data/{genome}.fasta",
        model_cfg="data/mammalia_softmasking_v2.yaml",
    output:
        gtf="results/tiberius/{genome}.gtf",
    params:
        seq_len=259992,
    container:
        "docker://quay.io/biocontainers/tiberius:1.1.6--pyhdfd78af_0"
    shell:
        # FIXME. python package doesn't get installed in biocontainer. Models
        # don't get shipped either.
        "python /usr/local/lib/python3.12/site-packages/tiberius/main.py "
        "--genome {input.fasta} "
        "--model_cfg {input.model_cfg} "
        "--out {output.gtf} "
        "--seq_len {params.seq_len} "
