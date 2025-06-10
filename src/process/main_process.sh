#!/bin/bash
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true

# source config file to get default parameters
CURRENT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CURRENT_DIR}/../config/config.sh"
source $sh_script

# Print help message
function print_help {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -s | --step <number>                    Step to execute (options: 1, 2, 3, 4, 5)"
    echo "                                          step 1. trim the adapter for the fastq files (R2 at default)"
    echo "                                          step 2. Mapping the trimmed fastq files to the reference genome (Use STAR)"
    echo "                                          step 3. Filter and sort the sam files"
    echo "                                          step 4. Generate the gene count matrix for all reads."
    echo "                                          step 5. Call SNPs and identify new reads from all reads"
    echo 
    echo "  -w | --working_dir <path>               Working directory to store all files"
    echo "                                          (Example: /data/qiyu/mRNAcap/data/241112_4SUtest_sample1and7)"
    echo "  -c | --cores <number>                   Number of CPU cores to use (default: 12)"
    echo "  -ids | --sample_ids <list>              Multiple sample IDs shoule be separated by comma withput sapce"
    echo "                                          (Example: [No1,S6hr1] or No1)"
    echo "  -fb | --filter_bam <true|false>         Whether to filter BAM files (default: true)"
    echo "  -rd | --remove_duplicates <true|false>  Whether to remove duplicate reads (default: false)"
    echo "  -rp | --rerun <true|false>              Rerun the pipeline with specific step (default: false)"
    echo "  --split_bam_type <type>                 Split Bam by READS or Chromosome groups (Provided: read | chromosome)"
    echo "  --fastq_suffix <suffix>                 FASTQ file suffix (default: .fastq.gz)"
    echo "  --seq_type <type>                       Sequencing type (Provided: Truseq | Nextera)"
    echo "  --species <species>                     Species (default: human)"
    echo "  --reference_folder <path>               Genome reference folder"
    echo "                                          (default: /broad/thechenlab/qiyu/mRNAcap/reference/Genome/GRCh38)"
    echo "  --genome_annotation <path>              Genome annotation GTF file"
    echo "                                          (default: ${reference_folder}/gencode.v39.primary_assembly.annotation.gtf)"
    echo "  --genome_fa <path>                      Genome FASTA file"
    echo "                                          (default: ${reference_folder}/GRCh38.primary_assembly.genome.fa)"
    echo "  -h | --help                             Display this help message and exit"
    echo 
}

# Parse argumentsr
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s | --step) step="$2"; shift 2;;
        -w | --working_dir) working_dir="$2"; shift 2;;
        -c | --cores) cores="$2"; shift 2;;
        -ids | --sample_ids) sample_ids="$2"; shift 2;;
        -fb | --filter_bam) filter_bam="$2"; shift 2;;
        -rd | --remove_duplicates) remove_duplicates="$2"; shift 2;;
        -rp | --rerun) rerun="$2"; shift 2;;
        --split_bam_type) rerun="$2"; shift 2;;
        --fastq_suffix) fastq_suffix="$2"; shift 2;;
        --seq_type) seq_type="$2"; shift 2;;
        --species) species="$2"; shift 2;;
        --reference_folder) reference_folder="$2"; shift 2;;
        --genome_annotation) genome_annotation="$2"; shift 2;;
        --genome_fa) genome_fa="$2"; shift 2;;
        -h|--help) print_help; exit 0;;
        *) echo "Unknown argument: $1"; print_help; exit 1;;
    esac
done

# check the parameters
if [ -z "$working_dir" ]; then
    echo "[ERROR] Working directory is not specified."
    exit 1
fi
if ! [[ "$step" =~ ^[1-5]$ ]]; then
    echo "[ERROR] --step must be one of: 1, 2, 3, 4, 5"
    exit 1
fi
if [[ "$split_bam_type" != "read" && "$split_bam_type" != "chromosome" ]]; then
    echo "[ERROR]  --split_bam_type must be 'read' or 'chromosome'"
    exit 1
fi
if [[ "$seq_type" != "Truseq" && "$seq_type" != "Nextera" ]]; then
    echo "[ERROR] --seq_type must be 'Truseq' or 'Nextera'"
    exit 1
fi
if [ ! -f "$genome_fa" ]; then
    echo "[ERROR] genome fasta not found at $genome_fa"
    exit 1
fi

# check the input fastq files for the sample ids
if [[ "$sample_ids" == *"["* && "$sample_ids" == *"]"* ]]; then
    sample_ids=${sample_ids//[\[\]]/}
fi
IFS=',' read -ra sample_ids <<< "$sample_ids"



##################### 1. trim the adapter for the fastq files (R2 at default) ######################
input_fastq_folder=${working_dir}/fastq
trim_input_folder=$input_fastq_folder
trim_output_folder="${working_dir}/1_trimmed"

if [ "$step" == 1 ]; then
    echo -e "\n[INFO] STEP 1: Start trimming the fastq files..."
    echo "[INFO] Sample IDs: ${sample_ids[@]}"

    echo input folder: $trim_input_folder
    echo output folder: $trim_output_folder

    cd $trim_input_folder

    fastq_files=()
    for sample_id in ${sample_ids[@]}; do
        files=($(find "$input_fastq_folder" -name "${sample_id}_*" | grep "_R2_.*${fastq_suffix}$"))
        if [ ${#files[@]} -ne 0 ]; then
            fastq_files+=("${files[@]}")
        fi
    done

    if [ ${#fastq_files[@]} -eq 0 ]; then
        echo "[ERROR] No fastq files found for all sample ids"
        exit 1
    fi

    if [ "$rerun" == "true" ]; then
        echo "[INFO] Re-run the trimming step"
        rm -rf $trim_output_folder
    fi

    mkdir -p $trim_output_folder
    chmod -R u+w $trim_output_folder

    # select the adapter based on the seq_type
    if [ "$seq_type" == "Nextera" ]; then
        adapter="${adapter:-CTGTCTCTTATACACATCT}"
    elif [ "$seq_type" == "Truseq" ]; then
        adapter="${adapter:-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT}"
    else
        echo "[ERROR] Invalid seq_type: $seq_type, should be Nextera or Truseq"
        exit 1
    fi
    echo "Selected adapter for $seq_type: $adapter"

    # trim the fastq files with cutadapt
    for file in "${fastq_files[@]}"; do
        echo -e "\nProcessing $file"
        file=$(basename "$file")
        output_file="${trim_output_folder}/${file}"
        tag_file="${output_file}.processing"

        if [ -f "$tag_file" ]; then
            rm -rf $output_file
            rm -rf $tag_file
        fi

        if [ -f "$output_file" ]; then
            echo "Skipping $file — already trimmed."
            continue
        fi

        touch "$tag_file"
        if $cutadapt -a "$adapter" -m 35 -j "$cores" \
            -o "$output_file" "${trim_input_folder}/${file}" --length 90; then
            echo "Trimming done for ${file}"
            rm -f "$tag_file"
        else
            echo "Trimming failed for ${file}"
            echo "Leaving tag file: $tag_file"
            exit 1
        fi
    done
fi



###################### 2. Mapping the trimmed fastq files to the reference genome (Use STAR) ######################
mapping_input_folder=${trim_output_folder}
STAR_output_folder="${working_dir}/2_alignment"

if [ "$step" == 2 ]; then
    echo -e "\n[INFO] STEP 2: Start alignment using STAR..."
    echo "[INFO] Sample IDs: ${sample_ids[@]}"

    cd $mapping_input_folder

    if [ "$rerun" == "true" ]; then
        echo "[INFO] Re-run the alignment step"
        rm -rf $STAR_output_folder
    fi

    mkdir -p $STAR_output_folder
    chmod -R u+w $STAR_output_folder

    if [ "$species" == "human" ]; then
        index="${reference_folder}/human_STAR_index/"
    elif [ "$species" == "mouse" ]; then
        index="${reference_folder}/mouse_STAR_index/"
    else
        echo "[ERROR] Invalid species: $species, should be human or mouse"
        exit 1
    fi

    # align read2 to the index file using STAR with default setting
    echo input folder: $mapping_input_folder
    echo output folder: $STAR_output_folder

    # start the alignment
    for file in "$mapping_input_folder"/*; do
        filename=$(basename "$file")

        for id in "${sample_ids[@]}"; do
            if [[ "$filename" == "${id}_"* ]]; then
                output_bam="${STAR_output_folder}/${filename}.Aligned.sortedByCoord.out.bam"
                progress_folder="${STAR_output_folder}/${filename}._STARtmp" 

                if [ -s "$output_bam" ] && [ ! -d "$progress_folder" ]; then
                    echo "Skipping $output_bam — output already exists and is non-empty."
                    break
                fi

                [ -f "$output_bam" ] && rm -f "$output_bam"
                [ -d "$progress_folder" ] && rm -rf "$progress_folder"

                echo -e "\nAligning $filename"
                $STAR --runThreadN "$cores" \
                    --outSAMtype BAM SortedByCoordinate \
                    --outSAMstrandField intronMotif \
                    --genomeDir "$index" \
                    --readFilesCommand zcat \
                    --readFilesIn "$file" \
                    --outFileNamePrefix "${STAR_output_folder}/${filename}."

                exit_code=$?
                output_bam="${STAR_output_folder}/${filename}.Aligned.sortedByCoord.out.bam"

                if [[ $exit_code -ne 0 ]]; then
                    echo "[ERROR] STAR failed for $filename (exit code $exit_code)"
                    exit 1
                elif [[ ! -s "$output_bam" ]]; then
                    echo "[ERROR] STAR finished but $output_bam is missing or empty"
                    exit 1
                else
                    echo "[INFO] Alignment done successfully for $filename"
                fi

                break
            fi
        done
    done

    # clean files
    rm -rf ${STAR_output_folder}/*Log.out
    rm -rf ${STAR_output_folder}/*Log.progress.out
fi



###################### 3. Filter and sort the sam files ######################
filtered_bam_folder="${working_dir}/3_bam"

if [ "$step" == 3 ]; then
    echo -e "\n[INFO] STEP 3: Start filtering bam files..."
    echo "[INFO] Sample IDs: ${sample_ids[@]}"

    if [ "$rerun" == "true" ]; then
        echo "[INFO] Re-run the filter and sort step"
        rm -rf $filtered_bam_folder
    fi

    mkdir -p $filtered_bam_folder
    chmod -R u+w $filtered_bam_folder

    echo input folder: $STAR_output_folder
    echo output folder: $filtered_bam_folder


    for file in $(ls "$STAR_output_folder"); do
        for id in "${sample_ids[@]}"; do
            if [[ "$file" == "${id}_"* && ( "$file" == *.bam || "$file" == *.sam ) ]]; then
                
                output_file="${filtered_bam_folder}/${file%.bam}_filtered.sorted.bam"
                output_file="${output_file%.sam}_filtered.sorted.bam"

                if [ -f "$output_file" ]; then
                    echo "Skipping $output_file — output already exists."
                    break
                fi

                echo -e "\nFiltering $file"
                $samtools view -bh -q 30 -F 4 "$STAR_output_folder/$file" | $samtools sort -@ 10 -o "$output_file"

                echo -e "Filtering done for $file"
                break
            fi
        done
    done

    # if remove_duplicates; then
    if [ "$remove_duplicates" == "true" ]; then
        echo -e "\nStart removing duplicates..."
        echo input folder: $filtered_bam_folder
        echo output folder: $filtered_bam_folder

        for file in $(ls "$filtered_bam_folder"); do
            for id in "${sample_ids[@]}"; do
                if [[ "$file" == "${id}_"* && "$file" == *_filtered.sorted.bam ]]; then
                    echo -e "\n[INFO] Removing duplicates for $file"

                    input_file="${filtered_bam_folder}/${file}"
                    output_file="${filtered_bam_folder}/${file%.bam}_rmdup.bam"

                    $samtools markdup -r "$input_file" "$output_file"
                    echo "Removing duplicates done for $file"

                    total_reads=$($samtools view -c "$input_file")
                    rmdup_reads=$($samtools view -c "$output_file")
                    echo "Total reads: $total_reads, reads after removing duplicates: $rmdup_reads"
                    
                    break
                fi
            done
        done

    fi

fi



###################### 4. Generate the gene count matrix ######################
counts_input_folder=${filtered_bam_folder}
counts_output_folder="${working_dir}/4_counts/all"

if [ "$step" == 4 ]; then
    echo -e "\n[INFO] STEP 4: Start generating gene count matrix..."
    echo "[INFO] Sample IDs: ${sample_ids[@]}"

    if [ "$remove_duplicates" == "true" ]; then
        echo "[INFO] Use the filtered and duplicates removed bam files for generating count matrix"
        input_files=($(ls $filtered_bam_folder | grep "_filtered.sorted.rmdup.bam$"))
    elif [ "$filter_bam" == "true" ]; then
        echo "[INFO] Use the filtered bam files for generating count matrix"
        input_files=($(ls $filtered_bam_folder | grep "_filtered.sorted.bam$"))
    else
        echo "[INFO] Use the original bam files for generating count matrix"
        counts_input_folder=$STAR_output_folder
        input_files=($(ls $STAR_output_folder | grep ".bam$"))
    fi

    if [ ${#input_files[@]} -eq 0 ]; then
        echo "[ERROR] No matching BAM files found."
        exit 0
    fi

    if [ "$rerun" == "true" ]; then
        echo "[INFO] Re-run the gene count matrix generation step"
        rm -rf $counts_output_folder
    fi
    
    mkdir -p $counts_output_folder
    chmod -R u+w $counts_output_folder

    echo input folder: $counts_input_folder
    echo output folder: $counts_output_folder

    for file in ${input_files[@]}; do
        echo -e "\nFeature counting for $file"

        sample_name=""
        for id in ${sample_ids[@]}; do
            if [[ $file =~ ^${id}_S[0-9]+_ ]]; then
                echo "Processing sample: $file"
                sample_name=$id
                break
            fi
        done
        if [[ -z $sample_name ]]; then
            continue 
        fi

        output_file_name="${counts_output_folder}/${sample_name}_matrix_counts.txt"
        if [ $cores -gt 64 ]; then
            cores=64
        fi
        if [ -f $output_file_name ]; then
            echo "[INFO] Skip feature counting for all reads"
        else
            $featureCounts -T $cores \
                -a $genome_annotation \
                -F GTF \
                --largestOverlap \
                -t exon \
                -g gene_name \
                -o $output_file_name \
                $counts_input_folder/${file}

            if [ $? -ne 0 ]; then
                echo "[ERROR] Feature counting failed for all reads"
                exit 1
            fi

            cat $output_file_name | cut -f1,7- | sed 1d > ${counts_output_folder}/${sample_name}_matrix_counts_genename.txt
            echo "Feature counting done for $file"
        fi
    done
fi



###################### 5. Call SNPs from all reads  ######################
snp_input_folder=$filtered_bam_folder
snp_output_folder="${working_dir}/5_SNP"
new_reads_output_folder="${working_dir}/4_counts/new"

if [ "$step" == 5 ]; then
    echo -e "\n[INFO] STEP 5: Start calling SNPs and identify new reads..."
    echo "[INFO] Sample IDs: ${sample_ids[@]}"

    if [ "$remove_duplicates" == "true" ]; then
        echo "[INFO] Use the filtered and duplicates removed bam for calling SNPs"
        snp_input_files=($(ls $filtered_bam_folder | grep "_filtered.sorted.rmdup.bam$"))
    elif [ "$filter_bam" == "true" ]; then
        echo "[INFO] Use the filtered bam for calling SNPs"
        snp_input_files=($(ls $filtered_bam_folder | grep "_filtered.sorted.bam$"))
    else
        echo "[INFO] Use the original bam for calling SNPs"
        snp_input_folder=$STAR_output_folder
        snp_input_files=($(ls $STAR_output_folder | grep ".bam$"))
    fi

    if [ "$rerun" == "true" ]; then
        echo "[INFO] Re-run the SNP calling step"
        rm -rf $snp_output_folder
        rm -rf $new_reads_output_folder
    fi

    mkdir -p $snp_output_folder
    mkdir -p $new_reads_output_folder
    chmod -R u+w $snp_output_folder
    chmod -R u+w $new_reads_output_folder

    echo input folder: $snp_input_folder
    echo output folder: $snp_output_folder

    for file in ${snp_input_files[@]}; do
        sample_name=""
        for id in ${sample_ids[@]}; do
            if [[ $file =~ ^${id}_S[0-9]+_ ]]; then
                sample_name=$id
                break
            fi
        done
        if [[ -z $sample_name ]]; then
            continue 
        fi
        echo -e "\n[INFO] Processing sample: $sample_name"

        # Pre-defined output folders
        # 5.1 mpileup and SNP calling
        output_mpileup="${snp_output_folder}/${sample_name}_output.mpileup"
        output_SNP_file="${snp_output_folder}/${sample_name}_SNP.vcf"
        # 5.2 BAM splitting
        bam_file="${snp_input_folder}/${file}"
        split_bam_folder="${snp_input_folder}/${sample_name}_split_bam"
        # 5.3 Identify new reads
        out_new_reads="${snp_output_folder}/${sample_name}_new_reads.txt"
        # 5.4 Subset bam files based on new reads
        subset_new_read_bam="${snp_input_folder}/${sample_name}_new_reads.bam"
        # 5.5 Generate count matrix for the new reads
        output_file_name="${new_reads_output_folder}/${sample_name}_matrix_counts.txt"

        # define the starting step based on the existing files
        starting_step=1
        if [ -f $output_file_name ]; then
            echo "[INFO] All steps for $sample_name completed."
            continue
        elif [ -f $subset_new_read_bam ]; then
            echo "[INFO] Skip subset bam files for $sample_name."
            starting_step=5
        elif [ -f $out_new_reads ]; then
            echo "[INFO] Skip new reads extraction for $sample_name."
            starting_step=4
        elif [[ ! -f "${split_bam_folder}/failed" && -d "$split_bam_folder" ]]; then
            echo "[INFO] Skip BAM splitting for $sample_name."
            starting_step=3
        elif [ -s $output_SNP_file ]; then
            echo "[INFO] Skip mpileup and SNP calling for $sample_name."
            starting_step=2
        else 
            starting_step=1
        fi
    
        # 5.1 Call SNPs for all reads
        if [ $starting_step -le 1 ]; then
            echo -e "\n[INFO] 5.1 Call SNPs for all reads"
            echo "samtools mpileup"
            $samtools mpileup -d 0 -BQ0 -R -f $genome_fa $snp_input_folder/${file} > $output_mpileup
            echo "varscan mpileup2snp"
            $java -jar $varscan mpileup2snp $output_mpileup --strand-filter 0 > $output_SNP_file

            if [ -f "$output_SNP_file" ]; then
                echo "SNP calling success: $output_SNP_file"
                rm -f "$output_mpileup"
            else
                echo "[ERROR] SNP output file not generated: $output_SNP_file"
                exit 1
            fi
        fi

        # 5.2 Split the bam file by READ_NAME groups
        if [ $starting_step -le 2 ]; then
            rm -rf "$split_bam_folder"
            mkdir -p "$split_bam_folder"

            if [ ! -f "${bam_file}.bai" ]; then
                samtools index "$bam_file" || { echo "[ERROR] Indexing failed."; exit 1; }
            fi

            # split bam file by READS group or Chromosome group
            if [[ "$split_bam_type" == "read" ]]; then
                echo -e "\n[INFO] 5.2 Split BAM by unique READ groups..."
                export -f process_read 
                export bam_file split_bam_folder sample_name

                # Extract unique READ_NAME and batch them in groups of bam_batch_size
                samtools view "$bam_file" | cut -f 1 | sort | uniq | \
                    split -l $bam_batch_size - "${split_bam_folder}/readname_batch_" \
                    || { echo "[ERROR] Failed to split by READ_NAME"; exit 1; }

                # find "${split_bam_folder}" -name "readname_batch_*" | parallel -j "$cores" process_read {} "$bam_file" "$split_bam_folder" "$sample_name"
                find "${split_bam_folder}" -name "readname_batch_*" | \
                    parallel -j "$cores" --results "${split_bam_folder}/parallel_output" \
                    'process_read {} "$bam_file" "$split_bam_folder" "$sample_name" || touch "$split_bam_folder/failed"'

            elif [[ "$split_bam_type" == "chromosome" ]]; then
                echo -e "\n[INFO] 5.2 Split BAM by unique Chromosome groups..."
                export -f process_chromosome
                export bam_file split_bam_folder sample_name

                # Get list of chromosomes, excluding '*' (unmapped reads)
                chromosomes=$(samtools idxstats "$bam_file" | cut -f 1 | grep -v '\*')
                echo "$chromosomes" | \
                    parallel -j "$cores" --results "${split_bam_folder}/parallel_output" \
                    'process_chromosome {} "$bam_file" "$split_bam_folder" "$sample_name" || touch "$split_bam_folder/failed"'
                
                # Unmapped reads
                if samtools idxstats "$bam_file" | grep -q '^\*'; then
                    samtools view -b "$bam_file" '*' > "${split_bam_folder}/${sample_name}_unmapped.bam"
                fi
            fi

            # Check for any failed batches
            if [ -f "${split_bam_folder}/failed" ]; then
                echo "[ERROR] Some batches failed during processing."
                echo "Detailed output can be found in ${split_bam_folder}/parallel_output"
                exit 1
            else
                echo "[INFO] All batches processed successfully."
            fi
        fi

        # 5.3 Identify new reads for each split bam file
        if [ $starting_step -le 3 ]; then
            echo -e "\n[INFO] 5.3 Identify T2C mutations"
            $python ${py_script} ${split_bam_folder} ${output_SNP_file} ${sample_name} ${genome_fa} ${snp_output_folder} ${cores} ${mut_rate_filter} ${end_dist}

            if [ $? -ne 0 ]; then
                echo "[ERROR] Failed to call mutations for $sample_name"
                exit 1
            fi
            if [ -f $out_new_reads ]; then
                echo "New reads extracted for $sample_name"
            else
                echo "[ERROR] Expected output file '$out_new_reads' not found for $sample_name"
                exit 1
            fi
        fi

        # 5.4 Subset bam files based on new reads
        if [ $starting_step -le 4 ]; then
            echo -e "\n[INFO] 5.4 Subset bam files based on new reads"
            subset_new_read_sam="${snp_input_folder}/${sample_name}_new_reads.sam"
            read_names_file="${snp_output_folder}/${sample_name}_read_names.txt"

            if [ -f $out_new_reads ]; then
                awk 'NR>1 {print $1}' $out_new_reads > $read_names_file
                samtools view -h $bam_file | awk 'BEGIN {while(getline < "'"$read_names_file"'") read_names[$1]=1} /^@/ || read_names[$1]' > $subset_new_read_sam
                samtools view -Sb $subset_new_read_sam > $subset_new_read_bam
                rm -rf $subset_new_read_sam
                rm -rf $read_names_file
            else
                echo "[ERROR] New reads file not found for $sample_name"
                exit 1
            fi
        fi

        # 5.5 Generate count matrix for the new reads
        if [ $starting_step -le 5 ]; then
            echo -e "\n[INFO] 5.5 Generate count matrix for the new reads"
            if [ ! -f $subset_new_read_bam ]; then
                echo "[ERROR] Subsetted bam file for new reads not found for $sample_name"
                exit 1
            fi
            if [ $cores -gt 64 ]; then
                cores=64
            fi

            $featureCounts -T $cores \
                -a $genome_annotation \
                -F GTF \
                --largestOverlap \
                -t exon \
                -g gene_name \
                -o $output_file_name \
                $subset_new_read_bam

            if [ $? -ne 0 ]; then
                echo "[ERROR] Feature counting failed for new reads"
                exit 1
            fi

            cat $output_file_name | cut -f1,7- | sed 1d > ${new_reads_output_folder}/${sample_name}_matrix_counts_genename.txt
            echo "[INFO] Feature counting done for new reads"
        fi

        # Clean intermediate files
        if [ -f "${split_bam_folder}/failed" ]; then
            echo "[ERROR] Some BAM batches failed during processing."
            echo "[HINT] Please rerun step 5.2 and all following steps."
            exit 1
        else
            rm -rf "${split_bam_folder}"
            echo "[INFO] All steps finished successfully."
        fi

    done

    # 5.6 Combine all and new counts into combined file
    $python ${combine_script} ${working_dir}
fi