
#!/bin/bash
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true

CURRENT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CURRENT_DIR}/../config/config.sh"

index_dir="${reference_folder}/${species}_STAR_index"  
mkdir -p $index_dir

$STAR --runThreadN $cores \
    --runMode genomeGenerate \
    --genomeDir $index_dir \
    --genomeFastaFiles $genome_fa \
    --sjdbGTFfile $genome_annotation \
    --sjdbOverhang 150

