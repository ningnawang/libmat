import os
import subprocess

def run_ftetwild_script(input_file, output_log):
    ftetwild_script = f"/home/junanita/ninwang/fTetWild/build/FloatTetwild_bin -i {input_file} -l 1 > {output_log}"
    # subprocess.run(ftetwild_script, shell=True)

def check_and_run_ftetwild(folder):
    output_log_folder = os.path.join(folder, 'ftetwild_output_log')
    os.makedirs(output_log_folder, exist_ok=True)

    for file_name in os.listdir(folder):
        file_path = os.path.join(folder, file_name)

        # Check if the item is a file and has the '.obj' extension
        if os.path.isfile(file_path) and file_name.endswith('.obj'):
            input_name = os.path.splitext(file_name)[0]
            output_log = os.path.join(output_log_folder, f"{input_name}_ftetwild_output.txt")

            print(f"Running fTetWild script for file: {file_name}")
            run_ftetwild_script(file_path, output_log)

def process_folders(base_folder):
    for folder_name in os.listdir(base_folder):
        folder_path = os.path.join(base_folder, folder_name)

        # Check if the item is a directory
        if os.path.isdir(folder_path):
            print(f"Checking folder: {folder_path}")
            check_and_run_ftetwild(folder_path)

if __name__ == "__main__":
    # Replace the base folder path with your actual directory
    # base_folder = "/home/junanita/ninwang/models/abc_0020_obj_v00"
    process_folders(base_folder)
