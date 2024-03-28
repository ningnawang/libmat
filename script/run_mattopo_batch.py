import os
import subprocess

def find_ma_file(msh_filename, output_ma_folder):
    # Extracting the common part of the filename
    common_part = msh_filename[:7]

    # Constructing the ma file name pattern
    ma_filename_pattern = f"mat_{common_part}"

    # Searching for ma files in the specified folder
    ma_files = [file for file in os.listdir(output_ma_folder) if file.startswith(ma_filename_pattern) and file.endswith(".ma")]

    if ma_files:
        # If ma file(s) exist, print a message and skip reading msh file
        print(f"Ma file(s) found for {msh_filename}. Skipping reading.")
        return True
    else:
        # If ma file does not exist
        return False

def run_bash_script(folder_path, file_name):
    input_file = os.path.join(folder_path, file_name)

    output_folder = os.path.join(folder_path, "mattopo_output_log")
    os.makedirs(output_folder, exist_ok=True)
    output_log = os.path.join(output_folder, f"{file_name}_output.txt")

    bash_script = f"compute-sanitizer --tool memcheck ../build/bin/VolumeVoronoiGPU {input_file} 100 > {output_log}"
    print(f"${bash_script}")

    with open(output_log, 'w') as log_file:
        try:
            subprocess.run(bash_script, shell=True, timeout=2400)  # (40 minutes)
        # except subprocess.TimeoutExpired:
        except:
            print(f"Timeout/Failed for file: {input_file}. Skipping.")

def main(folder_path, output_ma_folder):
    # Check if the folder exists
    if not os.path.exists(folder_path):
        print(f"Error: Folder '{folder_path}' does not exist.")
        return

    # Iterate over files in the folder
    for file_name in os.listdir(folder_path):
        # skip if not .msh file
        if not file_name.endswith('.msh'):
            continue

        # skip if .ma file exists
        is_ma_exist = find_ma_file(file_name, output_ma_folder) # boolean
        if is_ma_exist:
            continue

        file_path = os.path.join(folder_path, file_name)
        # Check if the item is a file (not a subfolder) and has a '.msh' extension
        if os.path.isfile(file_path) and file_name.endswith('.msh'):
            print(f"Running bash script for file: {file_name}")
            run_bash_script(folder_path, file_name)

if __name__ == "__main__":
    folder_path = '../../models/abc_test_2048_70/'
    # folder_path = '../data/xurui/'
    output_ma_folder = '../out/mat/'
    print(folder_path)
    main(folder_path, output_ma_folder)
