import pandas as pd
import os
import glob
import sys

def combine_count_files(working_dir):
    all_counts_folder = os.path.join(working_dir, "4_counts")
    output_file = os.path.join(all_counts_folder, "combined_all_new_counts.txt")
    count_dirs = [os.path.join(all_counts_folder, "new"), os.path.join(all_counts_folder, "all")]
    
    data_frames = []

    for count_dir in count_dirs:
        if not os.path.isdir(count_dir):
            print(f"Warning: {count_dir} does not exist, skipping.")
            continue

        file_list = glob.glob(os.path.join(count_dir, "*_matrix_counts_genename.txt"))
        for filepath in file_list:
            filename = os.path.basename(filepath)
            sample_name = filename.replace("_matrix_counts_genename.txt", "")
            suffix = "_n" if "new" in count_dir else "_a"
            sample_name += suffix

            try:
                print(f"Processing file: {filepath}")
                df = pd.read_csv(filepath, sep="\t", header=None, names=["Geneid", sample_name])
                df = df.dropna()

                df = df[pd.to_numeric(df[sample_name], errors='coerce').notnull()]
                df[sample_name] = df[sample_name].astype(int)

                df = df.set_index("Geneid")
                data_frames.append(df)
            except Exception as e:
                print(f"Failed to load {filepath}: {e}")

    if not data_frames:
        print("No valid count files found. Exiting.")
        return

    combined_df = pd.concat(data_frames, axis=1).fillna(0).astype(int)
    combined_df.index.name = "Geneid"
    combined_df.to_csv(output_file, sep="\t")

    print(f"Combined file written to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python combine_counts.py <working_dir>")
        sys.exit(1)

    working_dir = sys.argv[1]
    combine_count_files(working_dir)
