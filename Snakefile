configfile: "config.yaml"


import os
import os.path
import sys
import glob


# common parameters
primary_contigs_list = "resources/primary_contigs.list"
reference_fasta = "/cvmfs/softdrive.nl/projectmine_sw/resources/Build38/hs38DH/hs38DH.fa"
reference_index = "/cvmfs/softdrive.nl/projectmine_sw/resources/Build38/hs38DH/hs38DH.fa.fai"
reference_dict = "/cvmfs/softdrive.nl/projectmine_sw/resources/Build38/hs38DH/hs38DH.fa.dict"

# coverage inputs
preprocessed_intervals = "resources/preprocessed_intervals.interval_list"

# manta inputs
manta_region_bed = "resources/primary_contigs_plus_mito.bed.gz"

manta_region_bed_index = "resources/primary_contigs_plus_mito.bed.gz.tbi"

# PESR inputs
sd_locs_vcf = "/cvmfs/softdrive.nl/projectmine_sw/resources/Build38/hs38DH/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

# MELT inputs
melt_standard_vcf_header = "resources/melt_standard_vcf_header.txt"

melt_bed_file = config["melt_dir"] + "add_bed_files/Hg38/Hg38.genes.bed"

# whamg inputs
wham_include_list_bed_file = "resources/wham_whitelist.bed"

# module metrics parameters
primary_contigs_fai = "resources/contig.fai"


# checking_directory_structure
input_dir = config["input_dir"]
out_dir = config["out_dir"]
manta_dir = config["manta_dir"]
melt_dir = config["melt_dir"]

# checking sample names and input files can be found and read
sample_name = config["sample_name"]

if config["input_file_type"] in ["cram"]:
    input_files = expand(input_dir + "{sample}.cram", sample=sample_name)
    bam_dir = out_dir + "bam/"
    bam_file = bam_dir + "{sample}.bam"
    bam_index = bam_dir + "{sample}.bam.bai"
else:
    input_files = expand(input_dir + "{sample}.bam", sample=sample_name)
    bam_dir = input_dir
    bam_file =input_dir + "{sample}.bam"
    bam_index = input_dir + "{sample}.bam.bai"

# creating output_directories
filtered_bam_dir = out_dir + "filtered_bam_MELT/"
counts_dir = out_dir + "counts/"
pesr_dir = out_dir + "pesr/"
whamg_dir = out_dir + "whamg/"
multiple_metrics_dir = out_dir + "gatk_metrics/"
final_melt_dir = out_dir + "melt/"
scramble_dir = out_dir + "scramble/"
sample_metrics_dir = out_dir + "sample_metrics/"
module00cgvcf_dir = out_dir + "module00c_gvcf/"


# pseudo-rule to collect all target files

## ONLY IF INPUT FILE TYPE IS CRAM
if config["input_file_type"] in ["cram"]:

    rule bam:
        input:
            expand(bam_dir + "{sample}.bam", sample=sample_name)
            + expand(bam_dir + "{sample}.bam.bai", sample=sample_name),


rule counts:
    input:
        expand(counts_dir + "{sample}.counts.tsv.gz", sample=sample_name)
        + expand(counts_dir + "condensed_counts.{sample}.tsv.gz", sample=sample_name)
        + expand(counts_dir + "{sample}.interval_list", sample=sample_name),


rule pesr:
    input:
        expand(pesr_dir + "{sample}.pe.txt.gz", sample=sample_name)
        + expand(pesr_dir + "{sample}.pe.txt.gz.tbi", sample=sample_name)
        + expand(pesr_dir + "{sample}.sr.txt.gz", sample=sample_name)
        + expand(pesr_dir + "{sample}.sr.txt.gz.tbi", sample=sample_name)
        + expand(pesr_dir + "{sample}.sd.txt.gz", sample=sample_name)
        + expand(pesr_dir + "{sample}.sd.txt.gz.tbi", sample=sample_name),


rule multiplemetrics:
    input:
        expand(
            multiple_metrics_dir + "{sample}.alignment_summary_metrics",
            sample=sample_name,
        )
        + expand(
            multiple_metrics_dir + "{sample}.insert_size_metrics", sample=sample_name
        )
        + expand(
            multiple_metrics_dir + "{sample}.quality_distribution_metrics",
            sample=sample_name,
        ),


rule wgsmetrics:
    input:
        expand(multiple_metrics_dir + "{sample}_wgs_metrics.txt", sample=sample_name),


rule manta:
    input:
        expand(out_dir + "manta/{sample}.manta.vcf.gz", sample=sample_name)
        + expand(out_dir + "manta/{sample}.manta.vcf.gz.tbi", sample=sample_name),


rule melt:
    input:
        expand(final_melt_dir + "{sample}.melt.vcf.gz", sample=sample_name)
        + expand(final_melt_dir + "{sample}.melt.vcf.gz.tbi", sample=sample_name),


rule scramble:
    input:
        expand(scramble_dir + "{sample}.scramble.vcf.gz", sample=sample_name)
        + expand(scramble_dir + "{sample}.scramble.vcf.gz.tbi", sample=sample_name),


rule whamg:
    input:
        expand(whamg_dir + "{sample}.wham.vcf.gz", sample=sample_name)
        + expand(whamg_dir + "{sample}.wham.vcf.gz.tbi", sample=sample_name),


rule samplemetrics:
    input:
        expand(
            sample_metrics_dir + "standardised_vcf/{sample}.manta.std.vcf.gz",
            sample=sample_name,
        )
        + expand(
            sample_metrics_dir + "standardised_vcf/{sample}.melt.std.vcf.gz",
            sample=sample_name,
        )
        + expand(
            sample_metrics_dir + "standardised_vcf/{sample}.scramble.std.vcf.gz",
            sample=sample_name,
        )
        + expand(
            sample_metrics_dir + "standardised_vcf/{sample}.wham.std.vcf.gz",
            sample=sample_name,
        )
        + expand(
            sample_metrics_dir + "vcf_metrics/manta_{sample}.vcf.tsv",
            sample=sample_name,
        )
        + expand(
            sample_metrics_dir + "vcf_metrics/melt_{sample}.vcf.tsv",
            sample=sample_name,
        )
        + expand(
            sample_metrics_dir + "vcf_metrics/scramble_{sample}.vcf.tsv",
            sample=sample_name,
        )
        + expand(
            sample_metrics_dir + "vcf_metrics/wham_{sample}.vcf.tsv",
            sample=sample_name,
        )
        + expand(
            sample_metrics_dir + "SR_metrics/{sample}.sr-file.tsv", sample=sample_name
        )
        + expand(
            sample_metrics_dir + "PE_metrics/{sample}.pe-file.tsv", sample=sample_name
        )
        + expand(
            sample_metrics_dir + "counts_metrics/{sample}.raw-counts.tsv",
            sample=sample_name,
        ),


rule haplotype:
    input:
        expand(module00cgvcf_dir + "{sample}.g.vcf.gz", sample=sample_name)
        + expand(module00cgvcf_dir + "{sample}.g.vcf.gz.tbi", sample=sample_name),


# main rules that send the outputs to the target rules defined via the command line
if config["input_file_type"] in ["cram"]:

    rule CramToBam:
        input:
            reference_fasta,
            reference_index,
            cram_file=input_dir + "{sample}.cram",
        output:
            bam_file=bam_dir + "{sample}.bam",
            bam_index=bam_dir + "{sample}.bam.bai",
        benchmark:
            "benchmarks/CramToBam/{sample}.tsv"
        conda:
            "envs/samtools.yaml"
        threads: 4
        resources:
            mem_mb=2000,
        shell:
            """
            samtools view -b -h -@ {threads} -T {input[0]} --write-index -o {output.bam_file}##idx##{output.bam_index} {input.cram_file}
            """


rule CollectCounts:
    input:
        reference_fasta,
        reference_index,
        reference_dict,
        preprocessed_intervals,
        bam_file=bam_file,
        bam_index=bam_index,
    output:
        counts_file=counts_dir + "{sample}.counts.tsv.gz",
    benchmark:
        "benchmarks/CollectCounts/{sample}.tsv"
    conda:
        "envs/gatk.yaml"
    params:
        sample="{sample}",
        temp=counts_dir + "{sample}.counts.tsv",
    resources:
        mem_mb=12000,
    shell:
        """
        gatk --java-options "-Xmx10024m" CollectReadCounts --input {input.bam_file} --read-index {input.bam_index} --reference {input[0]} -L {input[3]} --format TSV --interval-merging-rule OVERLAPPING_ONLY -O {params.temp} --disable-read-filter MappingQualityReadFilter 
        sed -ri "s/@RG\tID:GATKCopyNumber\tSM:.+/@RG\tID:GATKCopyNumber\tSM:{params.sample}/g" {params.temp}
        bgzip {params.temp}
        """


rule CondenseReadCounts:
    input:
        counts=rules.CollectCounts.output.counts_file,
    output:
        condensed_counts=counts_dir + "condensed_counts.{sample}.tsv.gz",
    benchmark:
        "benchmarks/CondenseReadCounts/{sample}.tsv"
    params:
        sample="{sample}",
        temp_in_rd=counts_dir + "{sample}.in.rd.txt.gz",
        temp_out_rd=counts_dir + "{sample}.out.rd.txt.gz",
        temp_ref_dict=counts_dir + "{sample}.ref.dict",
    resources:
        mem_mb=12000,
    conda:
        "envs/gatk.yaml"
    shell:
        """
        zcat {input.counts} | grep '^@' | grep -v '@RG' > {params.temp_ref_dict}
        zcat {input.counts} | grep -v '^@' | sed -e 1d | awk 'BEGIN{{FS=OFS="\t";print "#Chr\tStart\tEnd\tNA21133"}}{{print $1,$2-1,$3,$4}}' | bgzip > {params.temp_in_rd}
        tabix -0 -s1 -b2 -e3 {params.temp_in_rd}
        gatk --java-options -Xmx2g CondenseDepthEvidence -F {params.temp_in_rd} -O {params.temp_out_rd} --sequence-dictionary {params.temp_ref_dict} --max-interval-size 2000 --min-interval-size 101
        cat {params.temp_ref_dict} <(zcat {params.temp_out_rd} | awk 'BEGIN{{FS=OFS="\t";print "@RG\tID:GATKCopyNumber\tSM:{params.sample}\\nCONTIG\tSTART\tEND\tCOUNT"}}{{if(NR>1)print $1,$2+1,$3,$4}}') | bgzip > {output.condensed_counts}
        rm {params.temp_ref_dict}
        rm {params.temp_in_rd}*
        rm {params.temp_out_rd}*
        """


rule CountsToIntervals:
    input:
        counts=rules.CollectCounts.output.counts_file,
    output:
        counts_intervals=counts_dir + "{sample}.interval_list",
    benchmark:
        "benchmarks/CountsToIntervals/{sample}.tsv"
    shell:
        """
        zgrep "^@" {input.counts} > {output.counts_intervals}
        zgrep -v "^@" {input.counts} | sed -e 1d | awk -F "\t" -v OFS="\t" '{{print $1,$2,$3,"+","."}}' >> {output.counts_intervals}
        """


rule PESRCollection:
    input:
        reference_fasta,
        reference_dict,
        sd_locs_vcf,
        preprocessed_intervals,
        primary_contigs_list,
        bam_file=bam_file,
    output:
        PE_file=pesr_dir + "{sample}.pe.txt.gz",
        PE_file_index=pesr_dir + "{sample}.pe.txt.gz.tbi",
        SR_file=pesr_dir + "{sample}.sr.txt.gz",
        SR_file_index=pesr_dir + "{sample}.sr.txt.gz.tbi",
        SD_file=pesr_dir + "{sample}.sd.txt.gz",
        SD_index=pesr_dir + "{sample}.sd.txt.gz.tbi",
    benchmark:
        "benchmarks/PESRCollection/{sample}.tsv"
    conda:
        "envs/gatk.yaml"
    threads: 1
    resources:
        mem_mb=10000,
    params:
        sample="{sample}",
    shell:
        """
        gatk --java-options "-Xmx3250m" CollectSVEvidence -I {input.bam_file} --sample-name {params.sample} \
        -F {input[2]} -SR {output.SR_file} -PE {output.PE_file} -SD {output.SD_file} \
        --site-depth-min-mapq 6 --site-depth-min-baseq 10 -R {input[0]} -L {input[4]} \
        --seconds-between-progress-updates 300
        """


rule runMultipleMetrics:
    input:
        reference_fasta,
        reference_index,
        bam_file=bam_file,
    output:
        alignment_file=multiple_metrics_dir + "{sample}.alignment_summary_metrics",
        insert_file=multiple_metrics_dir + "{sample}.insert_size_metrics",
        quality_file=multiple_metrics_dir + "{sample}.quality_distribution_metrics",
    benchmark:
        "benchmarks/runMultipleMetrics/{sample}.tsv"
    resources:
        mem_mb=10000,
    conda:
        "envs/gatk.yaml"
    params:
        sample="{sample}",
        metrics_base=out_dir + "gatk_metrics",
    shell:
        """
        gatk --java-options -Xmx3250m CollectMultipleMetrics -I {input.bam_file} -O {params.metrics_base}/{params.sample} -R {input[0]} --ASSUME_SORTED true --PROGRAM null --PROGRAM CollectAlignmentSummaryMetrics \
        --PROGRAM CollectInsertSizeMetrics --PROGRAM CollectSequencingArtifactMetrics --PROGRAM CollectGcBiasMetrics --PROGRAM QualityScoreDistribution --METRIC_ACCUMULATION_LEVEL null --METRIC_ACCUMULATION_LEVEL SAMPLE \
        --seconds-between-progress-updates 300
        """


rule runWGSMetrics:
    input:
        rules.runMultipleMetrics.output.alignment_file,
        reference_fasta,
        reference_index,
        bam_file=bam_file,
    output:
        wgs_metrics_file=multiple_metrics_dir + "{sample}_wgs_metrics.txt",
    benchmark:
        "benchmarks/runWGSMetrics/{sample}.tsv"
    resources:
        mem_mb=12000,
    conda:
        "envs/gatk.yaml"
    params:
        sample="{sample}",
        metrics_base=out_dir + "gatk_metrics",
    shell:
        """
        #extract read length as variable from multiple metrics
        readLength=$(cat {input[0]} | sed '8q;d' | awk '{{printf "%0.0f", $16}}')

        gatk --java-options -Xmx3250m CollectWgsMetrics --INPUT {input.bam_file} --VALIDATION_STRINGENCY SILENT --REFERENCE_SEQUENCE {input[1]} \
        --READ_LENGTH $readLength --INCLUDE_BQ_HISTOGRAM true --OUTPUT {params.metrics_base}/{params.sample}_wgs_metrics.txt --USE_FAST_ALGORITHM true \
        --VERBOSITY WARNING """


rule runMantastep1:
    input:
        reference_fasta,
        reference_index,
        manta_region_bed,
        manta_region_bed_index,
        bam_file=bam_file,
        bam_index=bam_index,
    output:
        diploid_vcf=out_dir + "manta/{sample}/results/variants/diploidSV.vcf.gz",
    benchmark:
        "benchmarks/runMantastep1/{sample}.tsv"
    conda:
        "envs/manta.yaml"
    shadow: "shallow"
    threads: 1
    resources:
        mem_mb=16000,
    params:
        memGb=16,
        sample="{sample}",
        manta_run_dir=out_dir + "manta/{sample}",
    shell:
        """
        rm -rf {params.manta_run_dir}
        {config[manta_dir]}bin/configManta.py --bam {input.bam_file} --referenceFasta {input[0]} --runDir {params.manta_run_dir} --callRegions {input[2]}
        {params.manta_run_dir}/runWorkflow.py --mode local --jobs {threads} --memGb {params.memGb} --quiet
        """


rule runMantastep2:
    input:
        reference_fasta,
        reference_index,
        rules.runMantastep1.output.diploid_vcf,
    output:
        manta_vcf=out_dir + "manta/{sample}.manta.vcf.gz",
        manta_index=out_dir + "manta/{sample}.manta.vcf.gz.tbi",
    benchmark:
        "benchmarks/runMantastep2/{sample}.tsv"
    conda:
        "envs/manta.yaml"
    threads: 1
    resources:
        mem_mb=16000,
    params:
        sample="{sample}",
        manta_run_dir=out_dir + "manta/{sample}",
    shell:
        """
        python2 {config[manta_dir]}libexec/convertInversion.py /$CONDA_PREFIX/bin/samtools {input[0]} {input[2]} | bcftools reheader -s <(echo "{params.sample}") > {params.manta_run_dir}/results/variants/diploidSV_inv.vcf
        bgzip -c {params.manta_run_dir}/results/variants/diploidSV_inv.vcf > {output.manta_vcf}
        tabix -p vcf {output.manta_vcf}
        """


rule runWhamg:
    input:
        reference_fasta,
        reference_index,
        wham_include_list_bed_file,
        bam_file=bam_file,
        bam_index=bam_index,
    output:
        whamg_vcf=whamg_dir + "{sample}.wham.vcf.gz",
        whamg_index=whamg_dir + "{sample}.wham.vcf.gz.tbi",
    benchmark:
        "benchmarks/runWhamg/{sample}.tsv"
    resources:
        mem_mb=32000,
    shadow: "shallow"        
    conda:
        "envs/whamg.yaml"
    params:
        sample="{sample}",
        whamg_dir_param=whamg_dir,
        whamg_c="chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY",
    shell:
        """
63963 Bus error
        awk 'BEGIN{{FS=OFS="\t"}}{{printf(\"%07d\\t%s\\n\",NR,$1\":\"$2\"-\"$3)}}' {input[2]} |\
          while read -r line interval; do
            vcfFile="$line.wham.vcf.gz"
            whamg -c "{params.whamg_c}" -x {threads} -a {input[0]} -f {input.bam_file} -r $interval | bgzip -c > $vcfFile
            bcftools index -t $vcfFile
          done

        ls -1 *.wham.vcf.gz > vcf.list
        bcftools concat -a -Ov -f vcf.list | sed -e 's/^#CHROM\t.*/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{params.sample}/' -e 's/;TAGS=[^;]*;/;TAGS={params.sample};/' | bgzip -c > {output.whamg_vcf}

        bcftools index -t {output.whamg_vcf}

        """


rule CreateHighCoverageIntervals:
    input:
        rules.CollectCounts.output.counts_file,
    output:
        high_coverage_intervals=final_melt_dir + "{sample}_highCountIntervals.bed",
    benchmark:
        "benchmarks/CreateHighCoverageIntervals/{sample}.tsv"
    resources:
        mem_mb=16000,
    params:
        threshold=1000,
    shell:
        """
        zgrep -v -E "^@|^CONTIG|^#" {input[0]} | \
          awk 'BEGIN{{FS=OFS="\t"}}{{if ($4 > {params.threshold}) {{print $1, $2, $3}}}}' > {output.high_coverage_intervals}
        """


rule CreateFilteredBAMMELT:
    input:
        rules.CreateHighCoverageIntervals.output.high_coverage_intervals,
        reference_fasta,
        reference_index,
        bam_file=bam_file,
    output:
        filtered_bam=filtered_bam_dir + "{sample}_filtered.bam",
    benchmark:
        "benchmarks/CreateFilteredBAMMELT/{sample}.tsv"
    conda:
        "envs/gatk.yaml"
    threads: 1
    resources:
        mem_mb=36000,
    params:
        interval_padding=100,
    shell:
        """
        gatk  --java-options "-Xmx3250m" PrintReads -XL {input[0]} --interval-exclusion-padding {params.interval_padding} -I {input[3]} -R {input[1]} --seconds-between-progress-updates 300 -O /dev/stdout | \
        samtools view -h - | awk 'BEGIN{{FS=OFS="\t"}}{{gsub(/[BDHVRYKMSW]/, "N", $10);print}}' | samtools view -b1 - > {output.filtered_bam}
        samtools index {output.filtered_bam}
        """


rule runMELT:
    input:
        rules.runWGSMetrics.output.wgs_metrics_file,
        rules.runMultipleMetrics.output.alignment_file,
        rules.runMultipleMetrics.output.insert_file,
        rules.CreateFilteredBAMMELT.output.filtered_bam,
        reference_fasta,
        reference_index,
        melt_bed_file,
    output:
        sva_vcf=final_melt_dir + "{sample}/SVA.final_comp.vcf",
        line1_vcf=final_melt_dir + "{sample}/LINE1.final_comp.vcf",
        alu_vcf=final_melt_dir + "{sample}/ALU.final_comp.vcf",
    benchmark:
        "benchmarks/runMELT/{sample}.tsv"
    threads: 1
    resources:
        mem_mb=36000,
    conda:
        "envs/MELT.yaml"
    params:
        sample="{sample}",
        melt_results_dir=final_melt_dir + "/{sample}",
    shell:
        """
        #get the coverage, read_length and insert_size from previous wgs and multiple metrics steps for more accurate MELT result
        meanCoverage=$(cat {input[0]} | sed '8q;d' | awk '{{printf $2}}')
        readLength=$(cat {input[1]} | sed '8q;d' | awk '{{printf "%0.0f", $16}}')
        insertSize=$(cat {input[2]} | sed '8q;d' | awk '{{printf "%0.0f", $6}}')

        #create transposon_reference_list
        ls {config[melt_dir]}/me_refs/Hg38/*zip | sed 's/\*//g' > {params.melt_results_dir}/transposon_reference_list
        
        #run MELTv2.2.2
        java -Xmx32000m -jar {config[melt_dir]}/MELT.jar Single -bamfile {input[3]} -h {input[4]} -c $meanCoverage -r $readLength -e $insertSize -d 40000000 -t {params.melt_results_dir}/transposon_reference_list -n {input[6]} -w {params.melt_results_dir}
        """


rule FixMELTOutput:
    input:
        rules.runMELT.output.sva_vcf,
        rules.runMELT.output.line1_vcf,
        rules.runMELT.output.alu_vcf,
        melt_standard_vcf_header,
        bam_file=bam_file,
    output:
        melt_vcf=final_melt_dir + "{sample}.melt.vcf.gz",
        melt_index=final_melt_dir + "{sample}.melt.vcf.gz.tbi",
    benchmark:
        "benchmarks/FixMELTOutput/{sample}.tsv"
    resources:
        mem_mb=4000,
    conda:
        "envs/MELT.yaml"
    params:
        sample="{sample}",
        melt_results_dir=final_melt_dir + "/{sample}",
        vcf_sort="resources/vcf-sort.pl",
        temp_header="/{sample}_temp_header.txt",
        temp_vcf="/{sample}.temp_melt.vcf.gz",
    shell:
        """
        cat {input[0]} | grep "^#" > "{params.melt_results_dir}/{params.sample}.header.txt"
        cat {input[0]} | grep -v "^#" > "{params.melt_results_dir}/{params.sample}.sva.vcf"
        cat {input[1]} | grep -v "^#" > "{params.melt_results_dir}/{params.sample}.line1.vcf"
        cat {input[2]} | grep -v "^#" > "{params.melt_results_dir}/{params.sample}.alu.vcf"
        cat {params.melt_results_dir}/{params.sample}.header.txt {params.melt_results_dir}/{params.sample}.sva.vcf {params.melt_results_dir}/{params.sample}.line1.vcf {params.melt_results_dir}/{params.sample}.alu.vcf | perl {params.vcf_sort} -c | bgzip -c > {params.melt_results_dir}/{params.sample}.melt_fix.vcf.gz
        tabix -p vcf {params.melt_results_dir}/{params.sample}.melt_fix.vcf.gz

        vcf_text=$(bgzip -cd {params.melt_results_dir}/{params.sample}.melt_fix.vcf.gz)
        grep '^#CHR' <<<"$vcf_text" | sed 's|{{basename({input.bam_file}, ".bam")}}|{params.sample}|g' > {params.melt_results_dir}/{params.temp_header}
        grep -v '^#' <<<"$vcf_text" | sed 's/No Difference/No_Difference/g' >> {params.melt_results_dir}/{params.temp_header}
        FILEDATE=$(grep -F 'fileDate=' <<<"$vcf_text")
        cat "{input[3]}" {params.melt_results_dir}/{params.temp_header} | sed "2i$FILEDATE" | bgzip -c > {params.melt_results_dir}/{params.sample}.melt_fix.vcf.gz
        mv {params.melt_results_dir}/{params.sample}.melt_fix.vcf.gz {params.melt_results_dir}/{params.temp_vcf}
        bcftools reheader -s <( echo "{params.sample}") {params.melt_results_dir}/{params.temp_vcf} > {output.melt_vcf}
        rm {params.melt_results_dir}/{params.temp_vcf}
        tabix -p vcf {output.melt_vcf}
        """


rule ScramblePart1:
    input:
        reference_fasta,
        primary_contigs_list,
        bam_file=bam_file,
    output:
        clusters_file=scramble_dir + "{sample}.scramble_clusters.tsv.gz",
    benchmark:
        "benchmarks/ScramblePart1/{sample}.tsv"
    singularity:
        "docker://us.gcr.io/broad-dsde-methods/markw/scramble:mw-scramble-99af4c50"
    threads: 1
    resources:
        mem_mb=2000,
    params:
        sample="{sample}",
    shell:
        """
        # Identify clusters of split reads
        while read region; do
          time /app/scramble-gatk-sv/cluster_identifier/src/build/cluster_identifier -l -r "$region" -t {input[0]} {input.bam_file} \
            | gzip >> {output.clusters_file}
        done < {input[1]}
        """


rule ScramblePart2:
    input:
        rules.ScramblePart1.output.clusters_file,
        reference_fasta,
    output:
        scramble_vcf=scramble_dir + "{sample}.scramble.vcf.gz",
        scramble_index=scramble_dir + "{sample}.scramble.vcf.gz.tbi",
    benchmark:
        "benchmarks/ScramblePart2/{sample}.tsv"
    singularity:
        "docker://us.gcr.io/broad-dsde-methods/markw/scramble:mw-scramble-99af4c50"
    threads: 8
    resources:
        mem_mb=16000,
    params:
        sample="{sample}",
        xDir=scramble_dir + "{sample}",
        clusterFile=scramble_dir + "{sample}/{sample}.scramble_clusters.tsv",
        ref="resources/blast_refs/",
    shell:
        """
        set -x
        mkdir -p {params.xDir}
        xDir=$PWD/{params.xDir}
        clusterFile=$PWD/{params.clusterFile}
        scrambleDir="/app/scramble-gatk-sv"
        meiRef=$scrambleDir/cluster_analysis/resources/MEI_consensus_seqs.fa

        gunzip -c {input[0]} > $clusterFile

        # Produce MEIs.txt
        Rscript --vanilla $scrambleDir/cluster_analysis/bin/SCRAMble.R --out-name $xDir/clusters --cluster-file $clusterFile --install-dir $scrambleDir/cluster_analysis/bin --mei-refs $meiRef --ref $PWD/{params.ref}ref --no-vcf --eval-meis --cores {threads}

        # create a header for the output vcf
        echo -e '##fileformat=VCFv4.3\n##reference={input[1]}\n##source=scramble' > $xDir/{params.sample}.tmp.vcf

        grep '^>' $meiRef | awk \
        '{{mei=toupper(substr($0,2)); if (mei=="L1") mei="LINE1"
          print "##ALT=<ID=INS:ME:" mei ",Description=\"" mei " element insertion\">"}}' >> $xDir/{params.sample}.tmp.vcf

        echo -e '##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Source algorithms\">\n##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Breakpoint strandedness [++,+-,-+,--]\">' >> $xDir/{params.sample}.tmp.vcf
        echo -e '##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate\">\n##INFO=<ID=MEI_START,Number=1,Type=Integer,Description=\"Start of alignment to canonical MEI sequence\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">' >> $xDir/{params.sample}.tmp.vcf

        blastdbcmd -db {params.ref}/ref -entry all -outfmt '##contig=<ID=%a,length=%l>' >> $xDir/{params.sample}.tmp.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{params.sample}" >> $xDir/{params.sample}.tmp.vcf

        # use awk to write the first part of an awk script that initializes an awk associative array
        # mapping MEI name onto its consensus sequence length
        awk \
        'BEGIN {{ FS=OFS="\t"; print "BEGIN{{" }}
         /^>/  {{if ( seq != "" ) print "seqLen[\"" seq "\"]=" len; seq = substr($0,2); len = 0}}
         !/^>/ {{len += length($0)}}
         END   {{if ( seq != "" ) print "seqLen[\"" seq "\"]=" len; print "}}"}}' $meiRef > $xDir/{params.sample}.awkScript.awk

        # write the rest of the awk script that transforms the contents of the *_MEIs.txt files into a VCF
        echo \
        'BEGIN{{ FS=OFS="\t" }}
        {{ if(FNR<2)next
          split($1,loc,":")
          start=loc[2]+1
          end=start+1
          len=seqLen[$2]-$10
          mei=toupper($2); if (mei=="L1") mei="LINE1"
          print loc[1],start,".","N","<INS:ME:" mei ">",int($6),"PASS",\
                "END=" end ";SVTYPE=INS;SVLEN=" len ";MEI_START=" $10 ";STRANDS=+-;CHR2=" loc[1] ";ALGORITHMS=scramble",\
                "GT","0/1" }}' >> $xDir/{params.sample}.awkScript.awk

        # transform the MEI descriptions into VCF lines
        awk -f $xDir/{params.sample}.awkScript.awk $xDir/clusters_MEIs.txt >> $xDir/{params.sample}.tmp.vcf

        #remove the tabs at the beginning of each line
        #sed -e 's/^\t*//' $xDir/{params.sample}.tmp.vcf > $xDir/{params.sample}.tmp.vcf

        # sort and index the output VCF
        bcftools sort -Oz <$xDir/{params.sample}.tmp.vcf >{output.scramble_vcf}
        bcftools index -ft {output.scramble_vcf}
        """


rule StandardiseVCF:
    input:
        rules.runMantastep2.output.manta_vcf,
        rules.FixMELTOutput.output.melt_vcf,
        rules.ScramblePart2.output.scramble_vcf,
        rules.runWhamg.output.whamg_vcf,
        primary_contigs_fai,
    output:
        standardised_vcf_manta=sample_metrics_dir
        + "standardised_vcf/{sample}.manta.std.vcf.gz",
        standardised_vcf_melt=sample_metrics_dir
        + "standardised_vcf/{sample}.melt.std.vcf.gz",
        standardised_vcf_scramble=sample_metrics_dir
        + "standardised_vcf/{sample}.scramble.std.vcf.gz",
        standardised_vcf_whamg=sample_metrics_dir
        + "standardised_vcf/{sample}.wham.std.vcf.gz",
    benchmark:
        "benchmarks/StandardiseVCF/{sample}.tsv"
    singularity:
        "docker://us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-01-24-v0.28.4-beta-9debd6d7"
    resources:
        mem_mb=16000,
    params:
        sample="{sample}",
        min_size=50,
        base_output_dir=sample_metrics_dir + "standardised_vcf/",
    shell:
        """
        svtk standardize --min-size {params.min_size} --contigs {input[4]} {input[0]} {params.base_output_dir}/{params.sample}.manta.std.vcf "manta"
        bgzip {params.base_output_dir}/{params.sample}.manta.std.vcf
        svtk standardize --min-size {params.min_size} --contigs {input[4]} {input[1]} {params.base_output_dir}/{params.sample}.melt.std.vcf "melt"
        bgzip {params.base_output_dir}/{params.sample}.melt.std.vcf
        svtk standardize --min-size {params.min_size} --contigs {input[4]} {input[2]} {params.base_output_dir}/{params.sample}.scramble.std.vcf "scramble"
        bgzip {params.base_output_dir}/{params.sample}.scramble.std.vcf
        svtk standardize --min-size {params.min_size} --contigs {input[4]} {input[3]} {params.base_output_dir}/{params.sample}.wham.std.vcf "wham"
        bgzip {params.base_output_dir}/{params.sample}.wham.std.vcf
        """


rule VCFMetrics:
    input:
        rules.StandardiseVCF.output.standardised_vcf_manta,
        rules.StandardiseVCF.output.standardised_vcf_melt,
        rules.StandardiseVCF.output.standardised_vcf_scramble,
        rules.StandardiseVCF.output.standardised_vcf_whamg,
        primary_contigs_list,
    output:
        vcf_metrics_manta=sample_metrics_dir + "vcf_metrics/manta_{sample}.vcf.tsv",
        vcf_metrics_melt=sample_metrics_dir + "vcf_metrics/melt_{sample}.vcf.tsv",
        vcf_metrics_scramble=sample_metrics_dir
        + "vcf_metrics/scramble_{sample}.vcf.tsv",
        vcf_metrics_whamg=sample_metrics_dir + "vcf_metrics/wham_{sample}.vcf.tsv",
    benchmark:
        "benchmarks/VCFMetrics/{sample}.tsv"
    singularity:
        "docker://us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-01-24-v0.28.4-beta-9debd6d7"
    resources:
        mem_mb=16000,
    params:
        sample="{sample}",
        types="DEL,DUP,INS,INV,BND",
        types_scramble="INS",
        base_samples_list=sample_metrics_dir + "vcf_metrics/",
    shell:
        """
        echo {params.sample} > {params.base_samples_list}/{params.sample}_name_list
        samplesList="{params.base_samples_list}/{params.sample}_name_list"
        svtest vcf {input[0]} {input[4]} $samplesList {params.types} "manta_{params.sample}" > {output.vcf_metrics_manta}
        svtest vcf {input[1]} {input[4]} $samplesList {params.types} "melt_{params.sample}" > {output.vcf_metrics_melt}
        svtest vcf {input[2]} {input[4]} $samplesList {params.types_scramble} "scramble_{params.sample}" > {output.vcf_metrics_scramble}
        svtest vcf {input[3]} {input[4]} $samplesList {params.types} "wham_{params.sample}" > {output.vcf_metrics_whamg}
        """


rule SRMetrics:
    input:
        rules.PESRCollection.output.SR_file,
    output:
        sr_metrics_file=sample_metrics_dir + "SR_metrics/{sample}.sr-file.tsv",
    singularity:
        "docker://us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-01-24-v0.28.4-beta-9debd6d7"
    resources:
        mem_mb=4000,
    params:
        sample="{sample}",
        base_samples_list=sample_metrics_dir + "SR_metrics/",
    shell:
        """
        echo {params.sample} > {params.base_samples_list}/{params.sample}_name_list
        samplesList="{params.base_samples_list}/{params.sample}_name_list"
        svtest sr-file {input[0]} $samplesList > {output.sr_metrics_file}
        """


rule PEMetrics:
    input:
        rules.PESRCollection.output.PE_file,
    output:
        pe_metrics_file=sample_metrics_dir + "PE_metrics/{sample}.pe-file.tsv",
    benchmark:
        "benchmarks/PEMetrics/{sample}.tsv"
    singularity:
        "docker://us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-01-24-v0.28.4-beta-9debd6d7"
    resources:
        mem_mb=4000,
    params:
        sample="{sample}",
        base_samples_list=sample_metrics_dir + "PE_metrics/",
    shell:
        """
        echo {params.sample} > {params.base_samples_list}/{params.sample}_name_list
        samplesList="{params.base_samples_list}/{params.sample}_name_list"
        svtest pe-file {input[0]} $samplesList > {output.pe_metrics_file}
        """


rule CountsMetrics:
    input:
        rules.CollectCounts.output.counts_file,
    output:
        counts_metrics_file=sample_metrics_dir
        + "counts_metrics/{sample}.raw-counts.tsv",
    benchmark:
        "benchmarks/CountsMetrics/{sample}.tsv"
    singularity:
        "docker://us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2024-01-24-v0.28.4-beta-9debd6d7"
    resources:
        mem_mb=4000,
    params:
        sample="{sample}",
    shell:
        """
        svtest raw-counts {input[0]} {params.sample} > {output.counts_metrics_file}
        """


rule Module00cGVCF:
    input:
        reference_fasta,
        bam_file=bam_file,
        interval_file=rules.CountsToIntervals.output.counts_intervals,
    output:
        GVCF=module00cgvcf_dir + "{sample}.g.vcf.gz",
        GVCF_index=module00cgvcf_dir + "{sample}.g.vcf.gz.tbi",
    benchmark:
        "benchmarks/Module00cGVCF/{sample}.tsv"
    conda:
        "envs/gatk.yaml"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx16G" HaplotypeCaller -R {input[0]} -I {input.bam_file} -L {input.interval_file} -O {output.GVCF} -ERC GVCF --native-pair-hmm-threads {threads} --interval-padding 100 
        """
