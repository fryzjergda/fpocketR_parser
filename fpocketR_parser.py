#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import glob
import os
import shutil
import subprocess
import zipfile
import pandas as pd

def argument_parser():
    """Parses command-line arguments for the script."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter, prog='fpocketR_parser')
    parser.add_argument("-f", "--file", required=True, dest="file", help="Input pdb file.")
    args = parser.parse_args()
    return args.file

def find_conda():
    """Finds the Conda executable path."""
    conda_path = os.environ.get('CONDA_EXE')
    if conda_path and os.path.isfile(conda_path):
        return conda_path

    conda_path = shutil.which('conda')
    if conda_path:
        return conda_path

    common_paths = [
        '/usr/local/miniconda3/bin/conda',
        '/usr/local/anaconda3/bin/conda',
        '/opt/conda/bin/conda',
        '/miniconda3/bin/conda',
        '/anaconda3/bin/conda'
    ]
    for path in common_paths:
        if os.path.isfile(path):
            return path
    return None

def run_fpocketR(inpdb):
    """Runs the fpocketR analysis on the specified pdb file and checks for the absence of pockets."""
    conda_executable = find_conda()
    if not conda_executable:
        print("Conda executable not found.")
        return False

    conda_env = 'fpocketR'
    fpocket_command = f'{conda_executable} run -n {conda_env} python -m fpocketR -pdb {inpdb} --ligand noll --offset 0 -o fpocket-R'
    
    try:
        result = subprocess.run(fpocket_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        print("fpocketR output:", result.stdout)
        print("fpocketR errors:", result.stderr)
        if "No pockets found" in result.stdout or "No pockets found" in result.stderr:
            print("No pockets were found by fpocketR.")
            return False  # Indicate that no pockets were found
        return True  # Indicate that the process was successful and pockets were likely found
    except subprocess.CalledProcessError as e:
        print("An error occurred while running fpocketR:", e)
        return False  # Indicate failure possibly due to execution error

def get_best_pockets(target_dir):
    """Identifies the best pockets from fpocketR analysis based on various metrics."""
    if not os.path.isdir(target_dir):
        print(f"Directory does not exist: {target_dir}")
        return

    os.chdir(target_dir)
    csv_files = glob.glob('*characteristics.csv')
    if not csv_files:
        print("No characteristics CSV file found.")
        return

    csv_file = csv_files[0]
    df = pd.read_csv(csv_file)
    relevant_data = df[['Pocket', 'Score', 'Drug score', 'SASA', 'Volume']]
    
    pockets = {
        'score': relevant_data.loc[relevant_data['Score'].idxmax()]['Pocket'],
        'drug_score': relevant_data.loc[relevant_data['Drug score'].idxmax()]['Pocket'],
        'sasa': relevant_data.loc[relevant_data['SASA'].idxmax()]['Pocket'],
        'volume': relevant_data.loc[relevant_data['Volume'].idxmax()]['Pocket']
    }
    
    os.chdir("./pockets")
    results = {}
    for metric, pocket in pockets.items():
        pocket_file = f"pocket{int(pocket)}_vert.pqr"
        if os.path.exists(pocket_file):
            df_pocket = parse_pqr(pocket_file)
            median_df = format_dataframe(df_pocket.median().to_frame().T)
            mean_df = format_dataframe(df_pocket.mean().to_frame().T)
            results[metric] = {'median': median_df, 'mean': mean_df}
            median_df.to_csv(f"{metric}_median.csv", index=False)
            mean_df.to_csv(f"{metric}_mean.csv", index=False)
    return results

def format_dataframe(df):
    """Formats the DataFrame by setting floating-point precision and ensuring integer pocket numbers."""
    df['pocket_number'] = df['pocket_number'].astype(int)
    for column in df.select_dtypes(include=['float64']).columns:
        df[column] = df[column].apply(lambda x: f"{x:.3f}")
    return df

def parse_pqr(file_path):
    """Parses PQR files to extract pocket numbers and coordinates."""
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                parts = line.split()
                x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
                pocket_number = int(parts[4])
                data.append([pocket_number, x, y, z])
    return pd.DataFrame(data, columns=['pocket_number', 'x', 'y', 'z'])

def move_and_compress_files(base_dir):
    """Moves CSV files from pockets to the current directory and compresses the fpocket-R directory."""
    print("Current path:", os.getcwd())
    source_dir_pattern = os.path.join(base_dir, '*_clean_out', 'pockets')
    source_dirs = glob.glob(source_dir_pattern)
    if not source_dirs:
        print(f"No directories found matching pattern: {source_dir_pattern}")
        return
    
    source_dir = source_dirs[0]
    csv_files = glob.glob(os.path.join(source_dir, '*.csv'))
    for file in csv_files:
        dest_file = os.path.join(os.getcwd(), os.path.basename(file))
        if os.path.exists(dest_file):
            os.remove(dest_file)  # Remove the file if it exists to allow overwriting
            print(f"Removed existing file at {dest_file}")
        shutil.move(file, dest_file)
        print(f"Moved {file} to {os.getcwd()}")

    zip_file = f"{base_dir}.zip"
    with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(base_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, os.path.join(base_dir, '..'))
                zipf.write(file_path, arcname)

    shutil.rmtree(base_dir)
    print(f"Removed directory {base_dir}")

if __name__ == '__main__':
    main_dir = os.getcwd()
    os.system("rm -r fpocket-R")
    infile = argument_parser()
    name = ''.join(infile.split(".")[:-1])
    

    if run_fpocketR(infile):
        char_dir = f"fpocket-R/{name}_clean_out"
        best_pockets_results = get_best_pockets(char_dir)
        print(best_pockets_results)

        base_directory = "./fpocket-R"
        os.chdir(main_dir)
        move_and_compress_files(base_directory)
    else:
        print("Processing terminated: No pockets found or fpocketR failed to execute properly.")
    os.system("rm -r fpocket-R")
