rule all:
    input:
        "plots_conda/quals.svg"

rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads_conda/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"



rule samtools_sort:
    input:
        "mapped_reads_conda/{sample}.bam"
    output:
        "sorted_reads_conda/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"



rule samtools_index:
    input:
        "sorted_reads_conda/{sample}.bam"
    output:
        "sorted_reads_conda/{sample}.bam.bai"
    conda:
      "envs/samtools.yaml"
    shell:
        "samtools index {input}"


"""
This is how you can tell snakemake to create a conda environment for the given rule, activate it
and deactivate it after the rule has finished.
Note that conda does not work with 'run', you must use either 'shell' (as seen here) or 'script'.

This is how to run the snakefile with conda:
    
    snakemake -s Snakefile_conda.smk --software-deployment-method conda --cores 10

IMPORTANT !!!
    The environemnt is not deleted after the snakemake script has finished.
    Find it via conda env list and delete (rm -rf <path>) it to save space.
"""



SAMPLES = ['A', 'B']

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads_conda/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads_conda/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls_conda/all.vcf"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"




rule plot_quals:
    input:
        "calls_conda/all.vcf"
    output:
        "plots_conda/quals.svg"
    script:
        "scripts/plot-quals.py"

