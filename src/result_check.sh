#!/bin/bash
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true

MAIN_DIR="/data/qiyu/mRNAcap/data"
FOLDER_NAME="250306_4SUdata_WalkUp"

input_fastq="${MAIN_DIR}/${FOLDER_NAME}/fastq"
output1="${MAIN_DIR}/${FOLDER_NAME}/1_trimmed"
output2="${MAIN_DIR}/${FOLDER_NAME}/2_alignment"
output3="${MAIN_DIR}/${FOLDER_NAME}/3_bam"
output4="${MAIN_DIR}/${FOLDER_NAME}/4_counts"
output5="${MAIN_DIR}/${FOLDER_NAME}/5_SNP"


# List all samples based on BAM files in output3
all_samples=($(ls "${output3}"/*sortedByCoord.out_filtered.sorted.bam 2>/dev/null | \
    xargs -n1 basename | cut -d '_' -f1 | sort -u))

echo -e "\n[INFO] Number of samples in total: ${#all_samples[@]}"


# 5.1. List samples with no SNP file (i.e., potentially empty or skipped)
missing_snp_files=()
echo -e "\n[INFO] Checking samples with missing SNP results..."
for sample in "${all_samples[@]}"; do
    snp_file="${output5}/${sample}_SNP.vcf"
    if [ ! -s "$snp_file" ]; then
        missing_snp_files+=("$sample")
    fi
done

if [ "${#missing_snp_files[@]}" -gt 0 ]; then
    echo "[WARN] SNP file missing for the following samples (step 5.1):"
    printf '%s\n' "${missing_snp_files[@]}"
else
    echo "[INFO] No samples failed with SNP file generating."
fi


# 5.2. Track samples with failed or empty BAM splitting
echo -e "\n[INFO] Checking samples with complete BAM batches..."

failed_bam_split_sample=()
empty_bam_split_sample=()

for split_bam_folder in "${output3}"/*_split_bam; do
    if [ -d "$split_bam_folder" ]; then
        sample_name=$(basename "$split_bam_folder" | sed 's/_split_bam$//')

        if [ -f "${split_bam_folder}/failed" ]; then
            failed_bam_split_sample+=("$sample_name")
        elif [ -z "$(find "$split_bam_folder" -name "*.bam" -type f 2>/dev/null)" ]; then
            # Folder exists but no .bam files
            empty_bam_split_sample+=("$sample_name")
        fi
    else
        # Folder does not exist
        sample_name=$(basename "$split_bam_folder" | sed 's/_split_bam$//')
        empty_bam_split_sample+=("$sample_name")
    fi
done

# Report missing or empty bam folders
if [ "${#empty_bam_split_sample[@]}" -gt 0 ]; then
    echo -e "\n[WARN] The following samples have missing or empty split BAM folders:"
    printf '%s\n' "${empty_bam_split_sample[@]}"
else
    echo "[INFO] All BAM split folders contain BAM files."
fi

# Report failures for bam split
if [ "${#failed_bam_split_sample[@]}" -gt 0 ]; then
    echo -e "\n[WARN] The following samples failed during BAM splitting (step 5.2):"
    printf '%s\n' "${failed_bam_split_sample[@]}"
else
    echo "[INFO] No samples failed during BAM splitting."
fi


# 5.3. Track samples with failed generation of new_reads.txt file
missing_new_reads_files=()
echo -e "\n[INFO] Checking samples with missing new_reads.txt..."
for sample in "${all_samples[@]}"; do
    txt_file="${output5}/${sample}_new_reads.txt"
    if [ ! -s "$txt_file" ]; then
        missing_new_reads_files+=("$sample")
    fi
done

if [ "${#missing_new_reads_files[@]}" -gt 0 ]; then
    echo "[WARN] new_reads.txt file missing for the following samples (step 5.3):"
    printf '%s\n' "${missing_new_reads_files[@]}"
else
    echo "[INFO] No samples failed with new_reads.txt file generating."
fi


# 5.4. Track samples with failed generation of new_reads.bam file
missing_new_bam_files=()
echo -e "\n[INFO] Checking samples with missing new_reads.bam..."
for sample in "${all_samples[@]}"; do
    bam_file="${output3}/${sample}_new_reads.bam"
    if [ ! -s "$bam_file" ]; then
        missing_new_bam_files+=("$sample")
    fi
done

if [ "${#missing_new_bam_files[@]}" -gt 0 ]; then
    echo "[WARN] new_reads.bam file missing for the following samples (step 5.4):"
    printf '%s\n' "${missing_new_bam_files[@]}"
else
    echo "[INFO] No samples failed with new_reads.bam file generating."
fi


# 5.5. Track samples with failed generation of new_reads_count.txt file
missing_new_count_files=()
echo -e "\n[INFO] Checking samples with missing new_reads_count.txt..."
for sample in "${all_samples[@]}"; do
    count_file="${output4}/new/${sample}_matrix_counts.txt"
    if [ ! -s "$count_file" ]; then
        missing_new_count_files+=("$sample")
    fi
done

if [ "${#missing_new_count_files[@]}" -gt 0 ]; then
    echo "[WARN] new_matrix_counts.txt file missing for the following samples (step 5.4):"
    printf '%s\n' "${missing_new_count_files[@]}"
else
    echo "[INFO] No samples failed with new_matrix_counts.txt file generating."
fi