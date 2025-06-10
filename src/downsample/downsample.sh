#!/bin/bash

sample=$1
folder=$2
downsample_rate=${2:-1} 
sample_file_type=${3:-"bam"} 
seed=${4:-100} 

if [ "$sample_file_type" == "bam" ]; then
    output_folder="/data/qiyu/mRNAcap/data/downsample/3_bam"
    mkdir -p "$output_folder"

    input_files=(/data/qiyu/mRNAcap/data/241220_4SUdata_22J5GLLT4/3_bam/${sample}_*.sorted.bam)

    if [[ ${#input_files[@]} -eq 0 ]]; then
        echo "Error: No files matching /data/qiyu/mRNAcap/data/241220_4SUdata_22J5GLLT4/3_bam/${sample}_*.sorted.bam found!"
        exit 1
    fi

    if [[ ${#input_files[@]} -gt 1 ]]; then
        echo "Warning: Multiple files found for ${sample}_*.sorted.bam. Using the first one."
    fi
    input_file="${input_files[0]}"
    echo "Using input file: $input_file"

    output_file="${output_folder}/${sample}_downsample_filtered.sorted.bam"
    seed_fraction="${seed}.${downsample_rate}"

    echo "Downsampling BAM file: $input_file with rate $downsample_rate..."
    samtools view -s "$seed_fraction" -b "$input_file" > "$output_file"

    if [[ -s "$output_file" ]]; then
        echo "Downsampled BAM file saved to $output_file"
    else
        echo "Error: Failed to generate downsampled BAM file."
        exit 1
    fi
elif [ $sample_file_type == "fastq" ]; then
    output_folder="/data/qiyu/mRNAcap/data/downsample/fastq"
    fastq_file="/data/qiyu/mRNAcap/data/241220_4SUdata_22J5GLLT4/fastq/${sample}_*R2_001.fastq.gz"
    output_file="${output_folder}/${sample}_downsampled.fastq"

    if [[ ! -f "$fastq_file" ]]; then
    echo "Error: Input FASTQ file $fastq_file does not exist!"
    exit 1
    fi

    # check if seqtk is installed
    if ! command -v seqtk &> /dev/null; then
    echo "Error: seqtk is not installed. Please install it first."
    exit 1
    fi

    # downsample FASTQ file
    echo "Downsampling $fastq_file at rate $downsample_rate..."
    seqtk sample -s100 "$fastq_file" "$downsample_rate" > "$output_file"

    # check if output file is generated
    if [[ -s "$output_file" ]]; then
    echo "Downsampled FASTQ file saved to $output_file"
    else
    echo "Error: Failed to generate downsampled FASTQ file."
    exit 1
    fi
fi