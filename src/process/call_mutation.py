import pysam
import pandas as pd
import os
import sys
import concurrent.futures
from tqdm import tqdm
import glob


def split_cigar_ops(cigartuples):
    """Split CIGAR operations into prefix, middle, and suffix"""
    prefix, suffix = [], []
    i, j = 0, len(cigartuples) - 1

    # Process prefix soft-clipping
    while i < len(cigartuples) and cigartuples[i][0] == 4:
        prefix.append(cigartuples[i])
        i += 1

    # Process suffix soft-clipping
    while j >= 0 and cigartuples[j][0] == 4:
        suffix.insert(0, cigartuples[j])
        j -= 1

    # Middle part
    middle = cigartuples[i:j+1]
    return prefix, middle, suffix


def fetch_reference_base(ref, chrom, pos, default='.'):
    """Fetch reference base from genome, with exception handling"""
    try:
        return ref.fetch(chrom, pos, pos + 1).upper()
    except Exception:
        return default


def generate_align_records(split_bam_path, ref_path, include_softclip=True, one_based_coords=True, include_deletions=False, show_progress=False):
    """
    Process BAM file to generate alignment records

    Parameters:
    - split_bam_path: Path to BAM file
    - ref_path: Path to reference genome
    - include_softclip: Whether to include soft-clipped parts
    - one_based_coords: Use 1-based (True) or 0-based (False) coordinate system
    - include_deletions: Whether to include deletions and skipped bases
    - show_progress: Whether to show progress bar

    Returns:
    - Generator yielding a dictionary per alignment record
    """
    bam = pysam.AlignmentFile(split_bam_path, "rb")
    ref = pysam.FastaFile(ref_path)
    cigar_dict = {i: j for i, j in zip(range(10), 'MIDNSHP=XB')}

    fetch_iterator = tqdm(bam.fetch()) if show_progress else bam.fetch(until_eof=True)

    for read in fetch_iterator:
        if read.is_unmapped:
            continue

        # Read basic info
        read_name, flag, chrom = read.query_name, read.flag, read.reference_name
        seq, qual = read.query_sequence, read.query_qualities or []
        ref_start, ref_end = read.reference_start, read.reference_end

        # Split CIGAR operations
        prefix, middle, suffix = split_cigar_ops(read.cigartuples)
        total_prefix_s = sum(op[1] for op in prefix)

        read_pos, ref_pos = 0, ref_start
        pos_adjust = 1 if one_based_coords else 0

        # Handle prefix soft-clipping
        if include_softclip:
            s_counter = 0
            for op, length in prefix:
                for i in range(length):
                    virtual_ref_pos = ref_start - total_prefix_s + s_counter
                    ref_base = fetch_reference_base(ref, chrom, virtual_ref_pos)

                    yield {
                        'READ_NAME': read_name,
                        'FLAG': flag,
                        'CHROM': chrom,
                        'READ_POS': read_pos,
                        'BASE': seq[read_pos] if read_pos < len(seq) else '.',
                        'QUAL': chr(qual[read_pos] + 33) if read_pos < len(qual) else '.',
                        'REF_POS': virtual_ref_pos + pos_adjust,
                        'REF': ref_base,
                        'OP': 'S'
                    }
                    read_pos += 1
                    s_counter += 1

        # Handle middle CIGAR operations
        for op, length in middle:
            op_chr = cigar_dict[op]
            for i in range(length):
                if op in (0, 7, 8):  # Match (M, =, X)
                    ref_base = fetch_reference_base(ref, chrom, ref_pos)

                    yield {
                        'READ_NAME': read_name,
                        'FLAG': flag,
                        'CHROM': chrom,
                        'READ_POS': read_pos,  
                        'BASE': seq[read_pos] if read_pos < len(seq) else '.',
                        'QUAL': chr(qual[read_pos] + 33) if read_pos < len(qual) else '.',
                        'REF_POS': ref_pos + pos_adjust,
                        'REF': ref_base,
                        'OP': op_chr
                    }
                    read_pos += 1
                    ref_pos += 1

                elif op == 1:  # Insertion (I)
                    ref_base = fetch_reference_base(ref, chrom, ref_pos)

                    yield {
                        'READ_NAME': read_name,
                        'FLAG': flag,
                        'CHROM': chrom,
                        'READ_POS': read_pos,
                        'BASE': seq[read_pos] if read_pos < len(seq) else '.',
                        'QUAL': chr(qual[read_pos] + 33) if read_pos < len(qual) else '.',
                        'REF_POS': ref_pos + pos_adjust,
                        'REF': ref_base,
                        'OP': 'I'
                    }
                    read_pos += 1

                elif op in (2, 3) and include_deletions:  # Deletion/Skipped (D/N)
                    ref_base = fetch_reference_base(ref, chrom, ref_pos)

                    yield {
                        'READ_NAME': read_name,
                        'FLAG': flag,
                        'CHROM': chrom,
                        'READ_POS': '.',
                        'BASE': '.',
                        'QUAL': '.',
                        'REF_POS': ref_pos + pos_adjust,
                        'REF': ref_base,
                        'OP': 'D' if op == 2 else 'N'
                    }
                    ref_pos += 1
                elif op in (2, 3) and not include_deletions:
                    ref_pos += 1

        # Handle suffix soft-clipping
        if include_softclip:
            s_counter = 0
            for op, length in suffix:
                for i in range(length):
                    virtual_ref_pos = ref_end + s_counter
                    ref_base = fetch_reference_base(ref, chrom, virtual_ref_pos)

                    yield {
                        'READ_NAME': read_name,
                        'FLAG': flag,
                        'CHROM': chrom,
                        'READ_POS': read_pos, 
                        'BASE': seq[read_pos] if read_pos < len(seq) else '.',
                        'QUAL': chr(qual[read_pos] + 33) if read_pos < len(qual) else '.',
                        'REF_POS': virtual_ref_pos + pos_adjust,
                        'REF': ref_base,
                        'OP': 'S'
                    }
                    read_pos += 1
                    s_counter += 1

    bam.close()
    ref.close()


def write_align_file(split_bam_path, ref_path, output_file, include_softclip=True, one_based_coords=True, include_deletions=True, show_progress=True):
    """
    Process BAM and write alignment records to file

    Parameters:
    - split_bam_path: BAM file path
    - ref_path: Reference genome path
    - output_file: Output file path
    - include_softclip: Whether to include soft-clipped regions
    - one_based_coords: Use 1-based coordinate system
    - include_deletions: Whether to include deletions/skipped bases
    - show_progress: Whether to display a progress bar
    """
    with open(output_file, "w") as fout:
        fout.write("READ_NAME\tFLAG\tCHROM\tREAD_POS\tBASE\tQUAL\tREF_POS\tREF\tOP\n")

        for record in generate_align_records(
            split_bam_path,
            ref_path,
            include_softclip=include_softclip,
            one_based_coords=one_based_coords,
            include_deletions=include_deletions,
            show_progress=show_progress
        ):
            fields = [
                record['READ_NAME'],
                str(record['FLAG']),
                record['CHROM'],
                str(record['READ_POS']),
                record['BASE'],
                record['QUAL'],
                str(record['REF_POS']),
                record['REF'],
                record['OP']
            ]
            fout.write("\t".join(fields) + "\n")

    print(f"Align file written to: {output_file}")


def call_mutation_from_records(SNP_df, record_iter, end_dist=0, mut_rate_filter=0.1):
    """
    Detect mutations from alignment records

    Parameters:
    - SNP_df: DataFrame of known SNPs
    - record_iter: Iterator of alignment records
    - end_dist: Distance from read ends to exclude
    - mut_rate_filter: Mutation ratio threshold

    Returns:
    - DataFrame of detected mutations, or None
    """
    print("Loading alignment records...")
    raw_df = pd.DataFrame.from_records(record_iter)
    
    # Convert numeric columns to numeric types
    for col in ['READ_POS', 'REF_POS', 'FLAG']:
        raw_df[col] = pd.to_numeric(raw_df[col], errors='coerce')
    
    # Filter out records with missing values in key columns
    raw_df.dropna(subset=['READ_NAME', 'CHROM', 'READ_POS', 'REF_POS', 'FLAG'], inplace=True)
    raw_df = raw_df.astype({'FLAG': int})  # Ensure FLAG is integer
    
    # Identify mutations (BASE != REF)
    df = raw_df[(raw_df['BASE'] != '.') & (raw_df['REF'] != '.') & (raw_df['BASE'] != raw_df['REF'])].copy()
    if df.empty:
        return None
    
    # Create chr_pos column for mutation identification
    df['REF_POS'] = df['REF_POS'].astype(int) + 1
    df['chr_pos'] = (
        df['CHROM'].astype(str) + '-' +
        df['REF_POS'].astype(int).astype(str) + '-' +
        df['REF'].astype(str) + '-' +
        df['BASE'].astype(str)
    )
    
    # Remove SNPs
    snp_set = set(SNP_df['chr_pos'].unique())
    df = df[~df['chr_pos'].isin(snp_set)]
    if df.empty:
        return None
    
    # Calculate end_point for each read
    end_pos_align = raw_df.groupby('READ_NAME')['READ_POS'].max().reset_index()
    end_pos_align.rename(columns={'READ_POS': 'end_point'}, inplace=True)
    
    # Join with end_point data
    df = pd.merge(df, end_pos_align, on='READ_NAME', how='left')
    
    # Filter by distance from read ends
    df = df[(df['READ_POS'] > end_dist) & (df['READ_POS'] < (df['end_point'] - end_dist))]
    if df.empty:
        return None
    
    # Calculate mutation statistics per read
    df['is_target'] = (
        ((df['FLAG'] == 0) & (df['REF'] == 'T') & (df['BASE'] == 'C')) |
        ((df['FLAG'] == 16) & (df['REF'] == 'A') & (df['BASE'] == 'G'))
    )
    
    # Group by read name and calculate mutation rate
    stats = df.groupby('READ_NAME').agg(
        mut_num=('READ_NAME', 'count'),
        target_mut_num=('is_target', 'sum')
    )
    stats['target_mut_ratio'] = stats['target_mut_num'] / stats['mut_num']
    
    # Filter reads by mutation rate
    good_reads = stats[stats['target_mut_ratio'] >= mut_rate_filter].index
    
    # Keep only target mutations from good reads
    filtered_df = df[df['READ_NAME'].isin(good_reads)]
    result = filtered_df[filtered_df['is_target']]
    
    if result.empty:
        return None
    
    # Select relevant columns
    result = result[['READ_NAME', 'FLAG', 'CHROM', 'READ_POS', 'BASE', 'QUAL', 
                    'REF_POS', 'REF', 'OP', 'chr_pos', 'end_point']]
    
    # Sort by READ_NAME and REF_POS as in the R code
    result = result.sort_values(['CHROM', 'REF_POS'])
    
    print("Done.")
    return result


def load_snp_data(snp_vcf_path):
    """
    Load SNP data from VCF file

    Parameters:
    - snp_vcf_path: Path to SNP VCF file

    Returns:
    - DataFrame of SNPs
    """
    try:
        SNP_df = pd.read_csv(snp_vcf_path, sep='\t')
    except:
        # Try with header=None if standard format fails
        try:
            SNP_df = pd.read_csv(snp_vcf_path, sep='\t', comment='#', header=None)
            SNP_df.columns = ['Chrom', 'Position', 'ID', 'Ref', 'Var', 'Qual', 'Filter', 'Info']
        except:
            raise ValueError(f"Could not parse SNP file: {snp_vcf_path}")

    if 'Chrom' in SNP_df.columns:
        SNP_df['chr_pos'] = (
            SNP_df['Chrom'].astype(str) + '-' +
            SNP_df['Position'].astype(str) + '-' +
            SNP_df['Ref'] + '-' +
            SNP_df['Var']
        )
    else:
        # Try to infer column names
        chrom_col = [col for col in SNP_df.columns if 'chrom' in col.lower() or 'chr' in col.lower()]
        pos_col = [col for col in SNP_df.columns if 'pos' in col.lower()]
        ref_col = [col for col in SNP_df.columns if 'ref' in col.lower()]
        alt_col = [col for col in SNP_df.columns if 'alt' in col.lower() or 'var' in col.lower()]

        if chrom_col and pos_col and ref_col and alt_col:
            chrom_col, pos_col = chrom_col[0], pos_col[0]
            ref_col, alt_col = ref_col[0], alt_col[0]

            # Ensure chromosome has 'chr' prefix
            if not SNP_df[chrom_col].astype(str).str.startswith('chr').all():
                SNP_df[chrom_col] = 'chr' + SNP_df[chrom_col].astype(str)

            SNP_df['chr_pos'] = (
                SNP_df[chrom_col].astype(str) + '-' +
                SNP_df[pos_col].astype(str) + '-' +
                SNP_df[ref_col] + '-' +
                SNP_df[alt_col]
            )

    return SNP_df


def process_bam_file(bam_path, SNP_df, ref_path, mut_rate_filter=0.1, end_dist=0):
    """
    Process a single BAM file and return the mutation results

    Parameters:
    - bam_path: Path to the BAM file
    - SNP_df: DataFrame of known SNPs
    - ref_path: Path to reference genome
    - mut_rate_filter: Mutation ratio threshold
    - end_dist: Distance from read ends to exclude

    Returns:
    - DataFrame of detected mutations, or None
    """
    try:
        print(f"Processing BAM file: {os.path.basename(bam_path)}")
        record_iter = generate_align_records(
            bam_path, 
            ref_path, 
            include_softclip=True, 
            one_based_coords=False, 
            include_deletions=False
        )
        
        result_df = call_mutation_from_records(
            SNP_df, 
            record_iter, 
            end_dist=end_dist, 
            mut_rate_filter=mut_rate_filter
        )
        
        return result_df
    except Exception as e:
        print(f"Error processing {bam_path}: {str(e)}")
        return None


def main():
    if len(sys.argv) < 4:
        print("Usage: python script.py <split_bam_path> <snp_vcf_path> <sample_name> <ref_path> [output_dir] [num_workers]")
        sys.exit(1)
    
    # Parse command line arguments
    split_bam_path = sys.argv[1]
    snp_vcf_path = sys.argv[2]
    sample_name = sys.argv[3]
    ref_path = sys.argv[4]
    output_dir = sys.argv[5]
    num_workers = int(sys.argv[6])
    mut_rate_filter = float(sys.argv[7])
    end_dist = int(sys.argv[8]) 

    print(f"Mut rate filter: {mut_rate_filter}")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load SNPs (shared among all processes)
    print(f"Loading SNP data from {snp_vcf_path}...")
    SNP_df = load_snp_data(snp_vcf_path)
    
    # Find all BAM files
    if os.path.isdir(split_bam_path):
        bam_files = glob.glob(os.path.join(split_bam_path, "*.bam"))
    else:
        bam_files = glob.glob(split_bam_path)
    
    if not bam_files:
        print(f"No BAM files found matching: {split_bam_path}")
        sys.exit(1)
    
    print(f"Found {len(bam_files)} BAM files to process")
    
    # Process BAM files in parallel
    all_results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Create a list of futures
        future_to_bam = {
            executor.submit(
                process_bam_file, 
                bam_file, 
                SNP_df, 
                ref_path,
                mut_rate_filter,
                end_dist
            ): bam_file for bam_file in bam_files
        }
        
        # Process completed futures as they complete
        for future in tqdm(concurrent.futures.as_completed(future_to_bam), total=len(bam_files)):
            bam_file = future_to_bam[future]
            try:
                result = future.result()
                if result is not None:
                    all_results.append(result)
                else:
                    print(f"No mutations detected in {os.path.basename(bam_file)}")
            except Exception as e:
                print(f"Error processing {bam_file}: {str(e)}")

    # Combine all results
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        combined_output_path = os.path.join(output_dir, f"{sample_name}_new_reads.txt")
        combined_df.to_csv(combined_output_path, sep='\t', index=False)
        print(f"Combined mutations output written to: {combined_output_path}")
        return combined_df
    else:
        print("No mutations detected in any of the BAM files.")
        return None

if __name__ == "__main__":
    main()
