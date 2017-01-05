
# Use BEDTools to get coverage of bam file
rule coverage:
    input: "{sample}.bam"
    output: "{sample}.coverage"
    shell: "{GenomeCoverageBed} -ibam {input} > {output}"
    
# Create index file for fast bam access
rule index:
    input: "{sample}.bam"
    output: "{sample}.bam.bai"
    log: "{sample}.bam.bai.sm.log"
    shell: "{SAMTOOLS} index {sample}"  

rule sort:
    input: "{sample}.bam"
    output: "{sample}.sorted.bam"
    log: "{sample}.sorted.log"
    shell: 
    "{SAMTOOLS} sort -T /tmp/{wildcards.sample}.sorted "
    " -O bam {output} {input}"  

rule samtobam:
    input: "{sample}.sam"
    output: "{sample}.bam"
    log: "{sample}.bam.sm.log"
    shell: """
    {SAMTOOLS} view -bS {sample} > {output}
    """






