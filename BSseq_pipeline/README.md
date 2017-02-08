# Snakemake workflow for BS-seq analysis

### How to run it
`snakemake --snakefile Snakemake.py`
### How to run it on cluster
e.g.
`snakemake --snakefile Snakemake.py --cluster "qsub -V -pe smp 8 -l h_rt=6:0:00 -l h_vmem=10G"`
### How to create a graph of rules
`snakemake --snakefile Snakemake.py -p -n --forceall --dag | dot -Tpdf > dag.pdf`

