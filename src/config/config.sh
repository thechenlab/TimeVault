#!/bin/bash

# Limit number of open files
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true

# ========== General Settings ==========
step=0
sample_ids=""
working_dir=""
cores="64"
filter_bam="true"
remove_duplicates="false"
rerun="false"
fastq_suffix=".fastq.gz"
seq_type="Nextera"
split_bam_type="read"
bam_batch_size=20000
mut_rate_filter=0.1
end_dist=0

# ========== Genome Reference ==========
species="human"
reference_folder="/data/qiyu/mRNAcap/reference"
genome_annotation="${reference_folder}/Genome/GRCh38/gencode.v39.primary_assembly.annotation.gtf"
genome_fa="${reference_folder}/Genome/GRCh38/GRCh38.primary_assembly.genome.fa"

# ========== Tools ==========
main_pkg="/mnt/share_pkgs"
python="/mnt/anaconda3/bin/python"
cutadapt="/mnt/anaconda3/bin/cutadapt"
java="/usr/bin/java"
samtools="/usr/local/bin/samtools"
STAR="/usr/bin/STAR"
featureCounts="/usr/bin/featureCounts"
varscan="${main_pkg}/VarScan.v2.3.9.jar"
picard="${main_pkg}/picard.jar"

# ========== Scripts ==========
main_script="/data/qiyu/mRNAcap/src"
py_script="${main_script}/process/call_mutation.py"
combine_script="${main_script}/process/combine_counts.py"
sh_script="${main_script}/process/batch_function.sh"