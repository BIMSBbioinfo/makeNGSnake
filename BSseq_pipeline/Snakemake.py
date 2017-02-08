#!/usr/bin/env python3.5
import os

# set config file
configfile: "./config.json"

# set working directory
WORKDIR = config["misc"]["work_dir"]


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
        [[expand(WORKDIR+"/Raw/fastq/fastqc/{sample}_fastqc.html",sample=config["unit"][sampleid]) for sampleid in config["samples"]]],
        [[expand(WORKDIR+"/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2_SE_report.html",sampleID=sampleid ,sample=config["unit"][sampleid]) for sampleid in config["samples"]]]
#        [[expand({{WORKDIR}}/Cleaned/Trimgalore/{sampleID}/fastqc/{sample}_trimmed_fastqc.html",sampleID=sampleid ,sample=config["unit"][sampleid]) for sampleid in config["samples"]]]


rule bismark_report:
    input:
        aln =  "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2_SE_report.txt",
        sp = "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated_splitting_report.txt",
        dd = "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.deduplication_report.txt",
        mbias = "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated.M-bias.txt",
        nuc = "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.nucleotide_stats.txt"
    output:
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2_SE_report.html"
    params:
        dir = "--dir {WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/"
    log: 
        "{WORKDIR}/Logs/bismark2report/{sample}.log"
    message: """--- Generate Bismark report."""
    shell:
        " {BISMARK2REPORT} {params} \
            --alignment_report {input.aln} \
            --splitting_report {input.sp} \
            --dedup_report {input.dd} \
            --mbias_report {input.mbias} \
            --nucleotide_report {input.nuc} \
            2> {log} "   


rule bismark_methylation_extractor:
    input:
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.deduplicated.bam"
    output:
        expand("{{WORKDIR}}/Mapped/Bismark/bowtie2/{{sampleID}}/non_directional/methylation_extracted/{type}_{strand}_{{sample}}_trimmed_bismark_bt2.deduplicated.txt.gz",type=["CHG","CHH","CpG"],strand=["OT","OB","CTOT","CTOB"]),
        expand("{{WORKDIR}}/Mapped/Bismark/bowtie2/{{sampleID}}/non_directional/methylation_extracted/{{sample}}_trimmed_bismark_bt2.deduplicated.{file}.gz",file=["bedGraph","bismark.cov","CpG_report.txt"]),
        
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated.M-bias.txt",
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated.M-bias_R1.png",
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted/{sample}_trimmed_bismark_bt2.deduplicated_splitting_report.txt"
    threads: 4
    params:
        se = "--single-end",
        gz = "--gzip",
        cReport = "--cytosine_report",
        bg = "--bedgraph",
        genomeFolder = "--genome_folder " + GENOME,
        outdir = "--output {WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/methylation_extracted"
    log: "{WORKDIR}/Logs/bismark_methylation_extraction/{sample}.log"
    message: """--- Extract Methylation Information."""
    shell:
        "{BISMARK_METHYLATION_EXTRACTOR} {params} --multicore {threads} {input} 2> {log}"


rule bismark_deduplication:
    input:
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.bam"
    output:
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.deduplicated.bam",
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.deduplication_report.txt"
    params:
        bam="--bam"
    log:
        "{WORKDIR}/Logs/deduplication/{sample}.log"
    message: """--- Deduplication step."""
    shell:
        "{DEDUPLICATE_BISMARK} {params} {input} 2> {log}"        
             
 

rule bismark_se_non_directional:
    input:
       "{WORKDIR}/Cleaned/Trimgalore/{sampleID}/{sample}_trimmed.fq.gz",
    output:
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.bam",
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2.nucleotide_stats.txt",
        "{WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/{sample}_trimmed_bismark_bt2_SE_report.txt"
    threads: 4
    params:
        N = "-N 1",
        L = "-L 20", 
        genomeFolder = "--genome_folder " + GENOME,
        outdir = "--output_dir {WORKDIR}/Mapped/Bismark/bowtie2/{sampleID}/non_directional/",
        nucCov = "--nucleotide_coverage",
        nonDir = "--non_directional ",
        pathToBowtie = "--path_to_bowtie "+ os.path.dirname(BOWTIE2) ,
        useBowtie2 = "--bowtie2 "
    log:
        "{WORKDIR}/Logs/bismark_mapping/bowtie2/{sample}.log"
    message: """--- Bismark bisulfite mapping to genome {VERSION}."""
    shell:
        "{BISMARK} {params} --multicore {threads} {input} 2> {log}"



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
        '{WORKDIR}/Logs/bismark_genome_preparation_'+VERSION+'.log'
    message: """--- Bisulfite conversion of genome: {VERSION}."""
    shell:
        "{BISMARK_GENOME_PREPARATION} {params} {input} 2> {log}"



rule fastqc_after_trimming:
    input:
        "{WORKDIR}/Cleaned/Trimgalore/{sampleID}/{sample}_trimmed.fq.gz"
    output:
    	"{WORKDIR}/Cleaned/Trimgalore/{sampleID}/fastqc/{sample}_trimmed_fastqc.html",
    	"{WORKDIR}/Cleaned/Trimgalore/{sampleID}/fastqc/{sample}_trimmed_fastqc.zip"
    params:
        outdir = "--outdir {WORKDIR}/Cleaned/Trimgalore/{sampleID}/fastqc/"
    log:
   	    "{WORKDIR}/Logs/fastqc_after_trimming/{sample}.log"
    message: """--- Quality check of trimmmed data with Fastqc."""
    shell:
      	"{FASTQC} {params.outdir} {input} 2> {log}" 

 

rule trimgalore_single_end:
    input:
        "{WORKDIR}/Raw/fastq/{sample}.fq"
    output:
        "{WORKDIR}/Cleaned/Trimgalore/{sampleID}/{sample}_trimmed.fq.gz",
        "{WORKDIR}/Cleaned/Trimgalore/{sampleID}/{sample}.fq_trimming_report.txt"
    params: 
        outdir = "--output_dir {WORKDIR}/Cleaned/Trimgalore/{sampleID}",
        phred = "--phred33",
        gz = "--gzip",
        cutadapt = "--path_to_cutadapt "+CUTADAPT,
        FivePrimeClip = "--clip_R1 10",
        ThreePrimeClip = "--three_prime_clip_R1 10"
    log:
        "{WORKDIR}/Logs/trimgalore/{sample}.trimgalore.log"
    message: 
        "--- Trimming of raw data with TrimGalore."
        "\n    Considering reads as single end."
    shell:
        "{TRIMGALORE} {params} {input} 2> {log}"



rule fastqc_raw:
    input:
	    "{WORKDIR}/Raw/fastq/{sample}.fq"
    output:
    	"{WORKDIR}/Raw/fastq/fastqc/{sample}_fastqc.html",
    	"{WORKDIR}/Raw/fastq/fastqc/{sample}_fastqc.zip"
    params:
        outdir = "--outdir {WORKDIR}/Raw/fastq/fastqc/"
    log:
   	    "{WORKDIR}/Logs/fastqc_raw/{sample}.log"
    message: """--- Quality check of raw data with Fastqc."""
    shell:
      	"{FASTQC} {params.outdir} {input} 2> {log}"
