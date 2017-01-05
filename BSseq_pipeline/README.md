# Snakemake workflow for BS-seq analysis

### How to run it
`snakemake --snakefile Snakemake.py`
### How to run it on cluster
e.g.
`snakemake --snakefile Snakemake.py --cluster "qsub -pe smp 10 -l h_rt=6:0:00 -l h_vmem=20G"`
