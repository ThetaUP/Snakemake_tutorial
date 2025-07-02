"""
All the software and the data is already there when you did the conda env.
"""

# info on this 'special rule are at the bottom of the script'
# FOR THE OTHER DEMO RULES TO WORK WITH SPECIFYING THE INPUTS IN THE CMD LINE YOU MUST
# COMMENT OUT THIS 'rule all'
rule all:
    input:
        "plots/quals.svg"

rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

"""
This is how you can run the rule above in order to generate multiple input files:

    snakemake -s Snakefile_1.smk --cores 1 mapped_reads/{A,B}.bam

Snakemake sees that you want to get to output files mapped_reads/A.bam and mapped_reads/B.bam.
Therefore it goes to the top of the pipeline and looks for two input files A.fastqc and B.fastqc
"""



rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

"""
See that in the rule above you can use the shell wildcars inside the snakemake scripts by using
'wildcars' (here 'wildcars.sample')
"""


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

"""
This is how you can visuallize the DAG of jobs (note that
the pipeline isn't actuallty executed here, we just do the dry run)

    snakemake -s Snakefile_1.smk sorted_reads/{A,B}.bam.bai --dag | dot -Tsvg > dag.svg
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


"""
In the rule above we use the expand() command to expand the filenames based on the list of possible files
(here ['A', 'B']).
This is a bit of a double edged sword, because when you hardcode this list inside the snakemake script there is
no way how you can then specify the samples via the command line. But there is a way to make this more flexible
via cofig files, but that will be in the 2d part of the tutorial.
"""


rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"


"""
In the rule above you can see how to outsorce the calculation to a script.

This is the (python) script:

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pysam import VariantFile

    quals = [record.qual for record in VariantFile(snakemake.input[0])]
    plt.hist(quals)

    plt.savefig(snakemake.output[0])


You can also run R (and some other languages) scripts.
"""


"""
There exists a special rule_all that must be placed as the first rule at the top. 
The idea behind this rule is that it specifies the output file you want to obtain,
so that you don't have to specify it in the command line.

Without the rule all we would run:
    snakemake -s Snakefile_1.smk -np plots/quals.svg

With the rule all we would run:
    snakemake -s Snakefile_1.smk -np
"""