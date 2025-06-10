# 🧬 TimeVault - mRNAcap Processing Pipeline
This repository contains scripts for processing 4SU-labeled mRNA capture data, including trimming, alignment, counting, SNP calling, and identifying new reads.
## 📁 Project Directory Structure
MAIN\_DIR/ (e.g. /data/qiyu/mRNAcap)
├── data/
│   └── 250306\_4SUdata\_WalkUp/
│       ├── \* submit\_jobs.sh        # Job submission script from `src` 
│       ├── 1\_trimmed/               # Trimmed FASTQ files
│       ├── 2\_alignment/             # STAR alignment results
│       ├── 3\_bam/                   # Sorted and indexed BAM files
│       ├── 4\_counts/                # Count matrices
│       ├── 5\_SNP/                   # SNP calling results
│       ├── fastq/                    # Input FASTQ files
│       └── logs/                     # Job logs
│
├── reference/
│   ├── Genome/
│   │   └── GRCh38/                   # Reference genome (FASTA, GTF)
│   └── human\_STAR\_index/           # STAR index files
│
└── src/
├── \* submit\_jobs.sh           
├── \* result\_check.sh          
├── annotation/                 # Gene annotation
├── config/                     # Configuration file
├── downsample/                 # Downsampling utilities
└── process/                    # Modular processing scripts
## 🚀 Steps 
1. **Clone the GitHub repository** to your local directory.
2. **Configure the package and reference directory** as outlined in `config.sh`.
3. **Copy `submit_jobs.sh`** into each relevant data folder.
4. **Run `submit_jobs.sh`** with your sample information, then execute:
   ```bash ./submit_jobs.sh```
5. **Locate the final combined count matrix** at:
   ```4_counts/combined_all_new_counts.txt```