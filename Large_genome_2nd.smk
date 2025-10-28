#!/usr/bin/env python3

import glob

# containers
tiberius = "docker://larsgabriel23/tiberius@sha256:796a9de5fdef73dd9148360dd22af0858c2fe8f0adc45ecaeda925ea4d4105d3"
bbmap = "docker://quay.io/biocontainers/bbmap:39.37--he5f24ec_0"  # new version for bp=t
# config
input_genomes = [
    "N_forsteri",
]


wildcard_constraints:
    genome="|".join(input_genomes),


# functions
def demux_files_for_genome(wildcards):
    pattern = f"results/{wildcards.genome}/partition/demux/genome.20.shred.*.fa"
    return sorted(glob.glob(pattern))


rule target:
    rule target:
    input:
        expand("results/{genome}/partition/demux/genome.20.shred.{chunk}.fa",
               genome=input_genomes,
               chunk=[0,1]) 
    shell:
        "echo 'done' > results/target_done.txt"


rule demuxbyname:
    input:
        "results/{genome}/partition/genome.20.shred.fa",
    output:
        directory("results/{genome}/partition/demux/"),
    log:
        "logs/partition/{genome}.demux.log",
    threads: 1
    resources:
        runtime=10,
        mem_mb=int(128e3),
    container:
        bbmap
    shell:
        "demuxbyname.sh -Xmx{resources.mem_mb}m "
        "namesubs=0 "
        "in={input} "
        "out={output}/genome.20.shred.%_.fa "
        "2>{log}"


rule shred:
    input:
        "results/{genome}/partition/genome.20.fa",
    output:
        "results/{genome}/partition/genome.20.shred.fa",
    log:
        "logs/partition/{genome}.shred.log",
    threads: 1
    resources:
        runtime=10,
        mem_mb=int(128e3),
    container:
        bbmap
    shell:
        "shred.sh -Xmx{resources.mem_mb}m "
        "length=500000000 "
        "overlap=1000000 "
        "equal=f "
        "in={input} "
        "out={output} 2>{log}"
