# 🧬 TimeVault
mRNAcap Processing Pipeline


This repository contains scripts for processing 4SU-labeled mRNA capture data, including trimming, alignment, counting, SNP calling, and identifying new reads.

## 📁 Project Structure
```
MAIN_DIR/ (e.g. /data/qiyu/mRNAcap)
├── data/
│   └── 250306_4SUdata_WalkUp/
│       ├── submit_jobs.sh          # Job submission script from `src` 
│       ├── 1_trimmed/              # Trimmed FASTQ files
│       ├── 2_alignment/            # STAR alignment results
│       ├── 3_bam/                  # Sorted and indexed BAM files
│       ├── 4_counts/               # Count matrices
│       ├── 5_SNP/                  # SNP calling results
│       ├── fastq/                  # Input FASTQ files
│       └── logs/                   # Job logs
│
├── reference/
│   ├── Genome/
│   │   └── GRCh38/                 # Reference genome (FASTA, GTF)
│   └── human_STAR_index/           # STAR index files
│
└── src/
│   ├── submit_jobs.sh           
│   ├── result_check.sh          
│   ├── annotation/                 # Gene annotation
│   ├── config/                     # Configuration file
│   ├── downsample/                 # Downsampling bam or fastq
│   └── process/                    # Modular processing scripts
```


## 🚀 Steps 
1. Clone the GitHub repository to your local directory.
2. Configure the package and reference directory as outlined in `config.sh`.
3. Copy `submit_jobs.sh` into run-specific folders, fill in your sample information, and execute it using: `bash ./submit_jobs.sh`.
4. Locate the final combined count matrix at: `4_counts/combined_all_new_counts.txt`.
