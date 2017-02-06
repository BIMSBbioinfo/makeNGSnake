#specify output dir
#OUTPUT_DIR="/home/agosdsc/test/bs_seq/Data/Raw"

# set config file
configfile: "./my_config2.json"

# define samples and unit
UNIT = config["unit"]
SAMPLES = config["samples"]

# set path for tools
FASTQC = config["tools"]["fastqc"]
TRIMGALORE = config["tools"]["trimgalore"]
CUTADAPT = config["tools"]["cutadapt"]
BOWTIE2 = config["tools"]["bowtie2"] 
BISMARK_GENOME_PREPARATION = config["tools"]["bismark_genome_preparation"]
BISMARK = config["tools"]["bismark"]
BISMARK_METHYLATION_EXTRACTOR = config["tools"]["bismark_methylation_extractor"]
DEDUPLICATE_BISMARK = config["tools"]["deduplicate_bismark"]
BISMARK2REPORT = config["tools"]["bismark2report"]


# define genome paths
GENOME = config["ref"]["genome"]
VERSION = config["ref"]["version"]
 
rule do_all:
    input:
        [[expand("Data/Raw/fastq/fastqc/{sample}_fastqc.html",sample=config["unit"][sampleid]) for sampleid in config["samples"]]],
        [[expand("Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2_SE_report.html",sampleID=sampleid ,sample=config["unit"][sampleid]) for sampleid in config["samples"]]]
#        [[expand("Data/Cleaned/Trimgalore/{sampleID}/fastqc/{sample}_trimmed_fastqc.html",sampleID=sampleid ,sample=config["unit"][sampleid]) for sampleid in config["samples"]]]


rule bismark_report:
    input:
        aln =  "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2_SE_report.txt",
        sp = "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated_splitting_report.txt",
        dd = "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.deduplication_report.txt",
        mbias = "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated.M-bias.txt",
        nuc = "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.nucleotide_stats.txt"
    output:
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2_SE_report.html"
    params:
        dir = "--dir Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/"
    log: 
        "Logs/bismark2report/{sample}.log"
    message: """--- Generate Bismark report."""
    shell:
        " {BISMARK2REPORT} {params} \
            --alignment_report {input.aln} \
            --splitting_report {input.sp} \
            --dedup_report {input.dd} \
            --mbias_report {input.mbias} \
            --nucleotide_report {input.nuc} \
            2> {log} "   

#rule extract_all:
#    input:
#       [[expand("Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated.bedGraph.gz",sampleID=sampleid ,sample=config["unit"][sampleid]) for sampleid in config["samples"]]]



rule bismark_methylation_extractor:
    input:
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.deduplicated.bam"
    output:
        expand("Data/Mapped/Bismark/bowtie2/{{sampleID}}/non_directional/methylation_extracted/{type}_{strand}_{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz",type=["CHG","CHH","CpG"],strand=["OT","OB","CTOT","CTOB"]),
        expand("Data/Mapped/Bismark/bowtie2/{{sampleID}}/non_directional/methylation_extracted/{{sample}}_trimmed_bismark_bt2.deduplicated.{file}.gz",file=["bedGraph","bismark.cov","CpG_report.txt"]),
        
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated.M-bias.txt",
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated.M-bias_R1.png",
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated_splitting_report.txt"
    threads: 4
    params:
        se = "--single-end",
        gz = "--gzip",
        cReport = "--cytosine_report",
        bg = "--bedgraph",
        genomeFolder = "--genome_folder " + GENOME,
        outdir = "--output Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted"
    log: "Logs/bismark_methylation_extraction/{sample}.log"
    message: """--- Extract Methylation Information."""
    shell:
        "{BISMARK_METHYLATION_EXTRACTOR} {params} --multicore {threads} {input} 2> {log}"

#rule bismark_all:
#    input:
#       [[expand("Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.bam",sampleID=sampleid ,sample=config["unit"][sampleid]) for sampleid in config["samples"]]]

rule bismark_deduplication:
    input:
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.bam"
    output:
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.deduplicated.bam",
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.deduplication_report.txt"
    params:
        bam="--bam"
    log:
        "Logs/deduplication/{sample}.log"
    message: """--- Deduplication step."""
    shell:
        "{DEDUPLICATE_BISMARK} {params} {input} 2> {log}"        
             
 

rule bismark_se_non_directional:
    input:
       "Data/Cleaned/Trimgalore/{sampleID}/{sample}_trimmed.fq.gz",
    output:
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.bam",
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.nucleotide_stats.txt",
        "Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2_SE_report.txt"
    threads: 4
    params:
#        multiCore = "--multicore {{threads}}",
        N = "-N 1",
        L = "-L 20", 
        genomeFolder = "--genome_folder " + GENOME,
        outdir = "--output_dir Data/Mapped/Bismark/bowtie2/{sampleID}/non_directional/",
        nucCov = "--nucleotide_coverage",
        nonDir = "--non_directional ",
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2 = "--bowtie2 "
    log:
        "Logs/bismark_mapping/bowtie2/{sample}.log"
    message: """--- Bismark bisulfite mapping to genome {VERSION}."""
    shell:
        "{BISMARK} {params} --multicore {threads} {input} 2> {log}"


import os

rule bismark_genome_preparation:
    input:
        "{GENOME}"
    output:
        "{GENOME}/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        "{GENOME}/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
    params:
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2 = "--bowtie2 ",
        verbose = "--verbose "
    log:
        'Logs/bismark_genome_preparation_'+VERSION+'.log'
    message: """--- Bisulfite conversion of genome: {VERSION}."""
    shell:
        "{BISMARK_GENOME_PREPARATION} {params} {input} 2> {log}"

#rule fastqc_trimmed_all:
#    input:
#        [[expand("Data/Cleaned/Trimgalore/{sampleID}/fastqc/{sample}_trimmed_fastqc.html",sampleID=sampleid ,sample=config["unit"][sampleid]) for sampleid in config["samples"]]]



rule fastqc_after_trimming:
    input:
        "Data/Cleaned/Trimgalore/{sampleID}/{sample}_trimmed.fq.gz"
    output:
    	"Data/Cleaned/Trimgalore/{sampleID}/fastqc/{sample}_trimmed_fastqc.html",
    	"Data/Cleaned/Trimgalore/{sampleID}/fastqc/{sample}_trimmed_fastqc.zip"
    params:
        outdir = "--outdir Data/Cleaned/Trimgalore/{sampleID}/fastqc/"
    log:
   	    "Logs/fastqc_after_trimming/{sample}.log"
    message: """--- Quality check of trimmmed data with Fastqc."""
    shell:
      	"{FASTQC} {params.outdir} {input} 2> {log}" 

 

rule trimgalore_single_end:
    input:
        "Data/Raw/fastq/{sample}.fq"
    output:
        "Data/Cleaned/Trimgalore/{sampleID}/{sample}_trimmed.fq.gz",
        "Data/Cleaned/Trimgalore/{sampleID}/{sample}.fq_trimming_report.txt"
    params: 
        outdir = "--output_dir Data/Cleaned/Trimgalore/{sampleID}",
        phred = "--phred33",
        gz = "--gzip",
        cutadapt = "--path_to_cutadapt "+CUTADAPT,
        FivePrimeClip = "--clip_R1 10",
        ThreePrimeClip = "--three_prime_clip_R1 10"
    log:
        "Logs/trimgalore/{sample}.trimgalore.log"
    message: 
        "--- Trimming of raw data with TrimGalore."
        "\n    Considering reads as single end."
    shell:
        "{TRIMGALORE} {params} {input} 2> {log}"

#rule fastqc_raw_all:
#    input:
#        [[expand("Data/Raw/fastq/fastqc/{sample}_fastqc.html",sample=config["unit"][sampleid]) for sampleid in config["samples"]]]


rule fastqc_raw:
    input:
	    "Data/Raw/fastq/{sample}.fq"
    output:
    	"Data/Raw/fastq/fastqc/{sample}_fastqc.html",
    	"Data/Raw/fastq/fastqc/{sample}_fastqc.zip"
    params:
        outdir = "--outdir Data/Raw/fastq/fastqc/"
    log:
   	    "Logs/fastqc_raw/{sample}.log"
    message: """--- Quality check of raw data with Fastqc."""
    shell:
      	"{FASTQC} {params.outdir} {input} 2> {log}"
