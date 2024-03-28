import os
import subprocess

def run_ftetwild_script(folder_path, file_name):
    input_file = os.path.join(folder_path, file_name)
    
    output_folder = os.path.join(folder_path, "ftetwild_output_log")
    os.makedirs(output_folder, exist_ok=True)
    output_log = os.path.join(output_folder, f"{file_name}_ftetwild_output.txt")

    ftetwild_script = f"/home/junanita/ninwang/fTetWild/build/FloatTetwild_bin -i {input_file} -l 0.5 > {output_log}"
    try:
        subprocess.run(ftetwild_script, shell=True, timeout=300)  # 300 seconds (5 minutes)
    except subprocess.TimeoutExpired:
        print(f"Timeout expired for file: {input_file}. Skipping.")

def main(folder_path):
    # Check if the folder exists
    if not os.path.exists(folder_path):
        print(f"Error: Folder '{folder_path}' does not exist.")
        return

    # Iterate over files in the folder
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)

        # Check if the '__sf.obj' file exists then continue
        if file_name.endswith('__sf.obj'):
            print(f"skip, {file_name} not target")
            continue

        msh_file_path = file_path + '_.msh'
        if not os.path.exists(msh_file_path) and file_name.endswith('.obj'):
            print(f"Running Tetwild script for file: {file_name}")
            run_ftetwild_script(folder_path, file_name)

if __name__ == "__main__":
    # Replace '../../models/abc_test_2048_100/' with the actual folder path containing your files
    # folder_path = '../../models/abc_test_2048_100/'
    folder_path = '../data/xurui/'
    main(folder_path)
