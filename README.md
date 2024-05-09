# gnomAD-SV_Module00
## **About**
This is the documentation for a Snakemake workflow which executes the GatherSampleEvidence module of the gatk-sv pipeline (available at https://github.com/broadinstitute/gatk-sv)

## Requirements
Snakemake installed via Conda in a standalone environment

## Directory structure 
-	/envs : these contain the conda environments for certain rules
-	/resources : the hg38 GATK resources and BLAST reference is in here
-	/tools : this directory contains the Manta and MELT executables with the hg38.genes.bed
-	/benchmarks : contains runtime metrics for each sample in each module
-	config.yaml : file containing the sample names, input/output dir 
-	Snakefile : main file containing the workflow code

## Downloads
This contains the download links for the reference (hg38) fasta and vcf files


**Usage:**
snakemake *COMMAND* -p --use-conda --use-singularity -j 100 "-B <IN_DIR>:<OUT_DIR>" --executor slurm --default-resources
