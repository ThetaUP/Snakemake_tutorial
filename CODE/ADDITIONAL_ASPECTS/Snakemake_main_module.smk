# this is how to include another snakemake file into your script
include: 'Snakemake_read_mapping_module.smk'  # it would be better if that would be the absolute path

"""
You just run this snakefile normally: snakemake -s Snakemake_main_module.smk --cores 10
"""

SAMPLES = ['A', 'B']

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"



rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"


