#!/usr/bin/env python3.5
import os

# path to the config file
configfile: "./config.json"

# output directory
OUTPUT_DIR = config["output_dir"]

# working directory
#workdir: OUTPUT_DIR # it gives error
os.chdir(OUTPUT_DIR)

TEMP_DIR = config["temp_dir"]

# paths to tools
FASTQC = config["tools"]["fastqc"]
BISMARK = config["tools"]["bismark"]
BOWTIE = config["tools"]["bowtie"]
GENOME = config["ref"]["genome"]
GENOME_BISULFITE_INDEX = config["ref"]["genome_index"]
CUTADAPT=config["tools"]["cutadapt"]
SAMTOOLS=config["tools"]["samtools"]
GenomeCoverageBed = config["tools"]["genomeCoverageBed"]
Deduplicate_Bismark = config["tools"]["deduplicate_bismark"]

# parameters
TRIM_GALONE_LENGTH = config["params"]["TRIM_GALONE_LENGTH"]
SINGLE_END = config["params"]["SINGLE_END"]
ADAPTER = config["params"]["ADAPTER"]
BBDUK_ADAPTER_ALL = config["params"]["BBDUK_ADAPTER"]
BBDUK_ADAPTER_TRUESEQ = config["params"]["BBDUK_ADAPTER_TRUESEQ"]

# output files
OUTPUT_FILES = [
expand("{sample}_fastqc.html", sample=config["samples"]),
expand("{sample}.cutadapt.fq", sample=config["samples"]),
expand("{sample}.cutadapt_fastqc.html", sample=config["samples"]),
expand("{sample}.cutadapt_bismark_bt2.bam", sample=config["samples"]),
expand("{sample}.cutadapt_bismark_bt2.deduplicated.bam", sample=config["samples"])
]


# rules
rule all:
  input: OUTPUT_FILES


rule deduplication:
  input: "{sample}.cutadapt_bismark_bt2.bam"
  output: "{sample}.cutadapt_bismark_bt2.deduplicated.bam"
  log: "{sample}.cutadapt_bismark_bt2.deduplication_report.txt"
  shell: "{Deduplicate_Bismark} --bam {input}"


rule bismark:
  input:
    genome = "/data/akalin/Base/Genomes/hg19/",
    file = "{sample}.cutadapt.fq"
  output:
    bam = "{sample}.cutadapt_bismark_bt2.bam"
  params:
    N = "1",
    L = "20",
  log:
    "{sample}_bismark_bt2_SE_report.txt"
  shell:
    """
    cmd="{BISMARK} --genome_folder {input.genome} --bowtie2 -o {OUTPUT_DIR} {input.file} -N {params.N} -L {params.L}";
    $cmd;
    mydate=`date +%Y-%m-%d:%H:%M:%S`;
    echo $mydate > bismark_command.log;
    echo $cmd >> bismark_command.log;
    """

#this takes very long time
rule bismark_genome_preparation:
  input:
    "{GENOME}"
  output:
    "{GENOME_BISULFITE_INDEX}/CT_conversion/genome_mfa.CT_conversion.fa",
    "{GENOME_BISULFITE_INDEX}/GA_conversion/genome_mfa.GA_conversion.fa"
  log:
    "{sample}_bismark_genome_preparation.log"
  shell:
    "{BISMARK_GENOME_PREPARATION} --path_to_bowtie {BOWTIE} --verbose {input} 2> {log}"


rule fastqc_after_trimming:
  input:
    "{sample}.cutadapt.fq"
  output:
    "{sample}.cutadapt_fastqc.html",
    "{sample}.cutadapt_fastqc.zip"
  log:
    "{sample}.cutadapt_fastqc.log"
  message: """--- Quality check of trimmed data with Fastqc."""
  shell:
    "{FASTQC} {input} 2> {log};"

rule cutadapt:
  input:
    file="{sample}.fastq"
  output:
    file="{sample}.cutadapt.fq"
  log:
    "{sample}.cutadapt.log"
  params:
    e="0.2",
    m="21",
    O="6",
    a=ADAPTER
  message: """--- Trimming with Cutadapt."""
  shell:
    """
    {CUTADAPT} -e {params.e} -a {params.a} -m {params.m} -O {params.O} {input.file} > {output.file} 2> {log}
    """

rule fastqc:
  input:
    "{sample}.fastq"
  output:
    "{sample}_fastqc.html",
    "{sample}_fastqc.zip"
  log:
    "{sample}_fastqc.log"
  message: """--- Quality check of raw data with Fastqc."""
  shell:
     "{FASTQC} {input} 2> {log}"

