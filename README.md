# makeNGSnake: the workflow for analysing the Next Generation Sequencing data

### How to run it:
```sh
$ snakemake -s Snakefile.py
```

### Content
  - Downloading fastq files from the ENA database based on their SRA ids
 
### Todo's
Add rules for:
  - Quality control 
  - Mapping downloaded data
  
### Requirements
Snakemake version 3.3
