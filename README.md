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
The download links for the reference (hg38) fasta and vcf files for running the pipeline

**Reference fasta and vcf** - https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false


      1.	Homo_sapiens_assembly38.fasta
      2.	Homo_sapiens_assembly38.dbsnp138.vcf, other files are already included in the resources folder

      
**Manta** - v1.6.0 - https://github.com/Illumina/manta/releases

**MELT** - v2.2.0 - https://melt.igs.umaryland.edu 

Downloads these files and store the reference files in the resources directory and the executables in the tools directory

## Changing the path of the input, output directories in the config.yaml file
      a.	input_dir : path to the input files (ends with /)
      b.	out_dir : path to the output files (ends with /)
      c.    Please include all the sample names in the batch as mentioned
      d.    Specify the Manta and MELT executables directory

## Changing the path of the resources in the snakefile
  Once the reference files and the resources folders is downloaded, keep them in the working directoty and mention the path in the Snakefile. 
      This includes,
      
**common parameters**

primary_contigs_list =  "resources/primary_contigs.list"
reference_fasta = "resources/Homo_sapiens_assembly38.fasta"
reference_index = "resources/Homo_sapiens_assembly38.fasta.fai"
reference_dict = "resources/Homo_sapiens_assembly38.dict"

**coverage inputs**

preprocessed_intervals = "resources/preprocessed_intervals.interval_list"

**manta inputs**

manta_region_bed = "resources/primary_contigs_plus_mito.bed.gz"
manta_region_bed_index = "resources/primary_contigs_plus_mito.bed.gz.tbi"

**PESR inputs** 

sd_locs_vcf = "resources/Homo_sapiens_assembly38.dbsnp138.vcf"

**MELT inputs**

melt_standard_vcf_header = "resources/melt_standard_vcf_header.txt"
melt_bed_file = config["melt_dir"] + "add_bed_files/Hg38/Hg38.genes.bed"

**whamg inputs**

wham_include_list_bed_file = "resources/wham_whitelist.bed"

**module metrics parameters** 

primary_contigs_fai = "resources/contig.fai"

Mention the PATH of the resources folder if not saved in the working directory

## Usage
Before running the pipeline, make sure snakemake environment is activated

snakemake *COMMAND* -p --use-conda --use-singularity -j 100 "-B <IN_DIR>:<OUT_DIR>" --executor slurm --default-resources
Mention the input and output directory for the bind options of singularity

## Submodule 
| Command | Description | Job allocation with 10 samples | Storage allocation with 10 samples | Memory allocation with 10 samples | CPU time with 10 samples
|--- | --- | --- | --- | --- | --- |
|Sample metrics| Runs the complete pipeline | -j 110 | 950 GB | 40 GB | 30 hours |





