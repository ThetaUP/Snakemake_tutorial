configfile: "config.yaml"  # talk about this below


rule all:
    input:
        "plots/quals.svg"


def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]



rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs

    output:
        "mapped_reads/{sample}.bam"         # comment this and uncomment the following when you are at the point of ...
        #temp("mapped_reads/{sample}.bam")  # ... talking about temporary files
    
    threads: 8

    log:
        "logs/bwa_mem/{sample}.log"

    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"

    shell:
        """
        (bwa mem -R '{params.rg}' -t {threads} {input} | 
        samtools view -Sb - > {output}) 2> {log}
        """

"""
In this rule we see how we can use parameters for our shell commands. This is done via the 'params' statement.
Here the parameters are the read groups (note that 'sample' stemms from the get_bwa_map_input_fastqs function).
We then use that parameter in the shell command via -R '{params.rg}' (because we names our parameter 'rg').
Note that the definition of the parameter uses plain python syntax. Look below at the bcftools_call rule for another
example of the usage of params.

You can also see that one of the inputs to this rule is a function that returns some value from the config file. Why is this
function needed? For example sometimes we don't know a priori the input for a given rule. For example, the input might change
given the output of some previous rule. To account for that we could inset an if-else statement into the function.

Moreover this is an example of a rule that has two inputs. These two inputs come from different rules.


This also shows how you can save the log files for each sample in the given process.
Note that the log files are not cleared even if the job fails.


Finally we can make a temprary file (via temp()) that gets deleted once the snakemake job has finished.
An example for this would be that we want to delete the .sam files once the .bam files have been created.
"""


"""
In the rule above you specify the exact number of threads the process will be run with.
But when you run the Snakfile from the command line then you MUST specify the number of cores (via --cores).
For example, if we specify 8 threads via threads: and then run snakemake with 10 cores (--cores 10),
then 8 cores will to the bwa_map rule and the remaining cores will get somehow redistributed amongst the
other rules. Offcourse if threads > cores, then the number of threads will be scaled down to match the number of cores.
If we want to use all the available cores on a given machine then we can use the '--cores all' option.
"""


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"              # comment this and uncomment the one below when you are at ...
        #protected("sorted_reads/{sample}.bam")  # ... the point of talking about this
    shell:
        """
        samtools sort -T sorted_reads/{wildcards.sample} 
        -O bam {input} > {output}
        """


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"



rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    params:
        prior_mutation_rate=config["prior_mutation_rate"]

    log:
        "logs/bcftools_call/all.log"
    shell:
        """
        bcftools mpileup -f {input.fa} {input.bam} | 
        bcftools call -P '{params.prior_mutation_rate}' -mv - > {output} 2> {log}
        """

"""
In the rule above you can see that we don't specify a list of samples like before (SAMPLES = ['A', 'B']).
Instead we specify the samples in a config.yaml (.json format also works) file that looks like this:

    samples:
    A: data/samples/A.fastq
    B: data/samples/B.fastq

Then, at the very top of this Snakefile we specify: configfile: "config.yaml".
And in this rule we reffer to the elements of the config file as config["samples"].


Moreover, we specify the parameter directly from the config file (prior_mutation_rate=config["prior_mutation_rate"]) and then
we use that parameter in bcftools call (-P '{params.prior_mutation_rate}').


Finally, here we also save the log files. But note that no wildcard ({sample}) is used here for logging, 
since we only get one output (once vcf for the entire script).

Additionally we can say that we want an output file to be protected (protected()) from accidental removals.
"""


rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"