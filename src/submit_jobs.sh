#!/bin/bash
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true

# Submit sample information
MAIN_DIR="/data/qiyu/mRNAcap"
FOLDER_NAME='250504_4SUdata_WalkUp'
sample_list=(
    S121
    S122
    S123
    S241
    S242
    S243
    S481
    S482
    S483
    S721
    S722
    S723
)

# Requested resource and running parameters
SERVER="LOCAL" # "UGER" or "LOCAL"
CORE=32          
MEM=18G 
TIME=1080  
STEP=1 # 0-all; 1, 2, 3, 4, 5
RERUN="false"
MUT_RATE_FILTER=0.1


# Run processing
export STEP CORE RERUN SERVER MUT_RATE_FILTER sample_list
source "${MAIN_DIR}/src/process/submit_function.sh"
submit_function "$MAIN_DIR" "$FOLDER_NAME"