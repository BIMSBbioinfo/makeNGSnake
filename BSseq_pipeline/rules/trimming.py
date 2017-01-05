


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
    
rule bbduk2:
  input:
    file="{sample}.fastq"
  output:
    file="{sample}.bbduk.fq"
  params:
    a=BBDUK_ADAPTER,
    minlength="21",
    qtrim="r",
    trimq="20",
    ktrim="r",
    k="25",
    mink="11",
    hdist=1,
    overwrite="false"
  log:
    "{sample}.bbduk.log"
  message: """--- Trimming with Bbduk2."""
  shell:
    "{BBDUK} -Xmx1g in={input.file} out={output.file} minlength={params.minlength} qtrim={params.r} \
    trimq={params.trimq} ktrim={params.ktrim} k={params.k} mink={params.mink} ref={params.a} hdist={params.hdist} \
    overwrite={params.overwrite} 2> {log}"

rule trim_galone:
  input:
    "{sample}.fastq"
  output:
    "{sample}_trimmed.fq"
  params:
    length = TRIM_GALONE_LENGTH
  log:
    "{sample}.fastq_trimming_report.txt"
  message: """--- Trimming with Trim Galone!."""
  shell:
    "{TRIM_GALONE} {input} --output_dir {OUTPUT_DIR} --length {params.length} --path_to_cutadapt {CUTADAPT}"
    
