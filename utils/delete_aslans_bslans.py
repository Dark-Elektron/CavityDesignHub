import os
import glob

# Define the base directory
base_dir = r"D:\Dropbox\CavityDesignHub\PhD_Thesis\SimulationData\SLANS"

# Define the filenames to delete
files_to_delete = ["aslans.mtx", "bslans.mtx"]

# Walk through the directory
for root, dirs, files in os.walk(base_dir):
    for filename in files_to_delete:
        file_path = os.path.join(root, filename)
        if os.path.isfile(file_path):
            # print(file_path)
            try:
                os.remove(file_path)
                # print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Failed to delete {file_path}: {e}")

print("Completed file deletion.")