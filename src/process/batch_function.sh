#!/bin/bash


# 1. Split BAM file by Chromosome groups
process_chromosome() {
    local chromosome=$1
    local bam_file=$2
    local split_bam_folder=$3
    local sample_name=$4
    
    output_bam="${split_bam_folder}/${sample_name}_${chromosome}.bam"
    echo "Processing chromosome: $chromosome"
    
    samtools view -b "$bam_file" "$chromosome" > "$output_bam"

    # Validate the generated BAM file
    if [ ! -f "$output_bam" ]; then
        echo "Error: BAM file not created for chromosome $chromosome"
        return 1
    fi

    samtools quickcheck "$output_bam"
    if [ $? -ne 0 ]; then
        echo "Error: Invalid BAM file generated for chromosome $chromosome"
        return 1
    fi

    samtools index "$output_bam"

    return 0
}


# 2. Split BAM file by READS groups
process_read() {
    local batch_file=$1
    local bam_file=$2
    local split_bam_folder=$3
    local sample_name=$4

    batch_name=$(basename "$batch_file")
    output_bam="${split_bam_folder}/${sample_name}_${batch_name}.bam"
    echo "Processing batch: $batch_name"

    {
        samtools view -H "$bam_file" 
        samtools view "$bam_file" | grep -F -f "$batch_file" 
    } | samtools view -b -o "$output_bam" -

    # Check if BAM file was created and is valid
    if [ ! -f "$output_bam" ]; then
        echo "Error: BAM file not created for batch $batch_name"
        return 1
    fi

    samtools quickcheck "$output_bam"
    if [ $? -ne 0 ]; then
        echo "Error: Invalid BAM file generated for batch $batch_name"
        samtools fillmd -u "$output_bam" "$bam_file" > "${output_bam}.tmp"
        if [ $? -eq 0 ]; then
            mv "${output_bam}.tmp" "$output_bam"
            echo "Attempted to repair BAM file for batch $batch_name"
        else
            return 1
        fi
    fi

    return 0
}
