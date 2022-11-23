# convert matlab in folders to python

#"python matlab2python/matlab2python.py SSC_single_model.m -o SSC_single_model.py"
import os
import subprocess

dirr = fr'{os.getcwd()}\multipacting_codes'
print(dirr)

if os.path.exists(dirr):
    print("yeah")

for path, subdirs, files in os.walk(dirr):
    for name in files:
        if name.split(".")[-1] == 'm':
            file_dirname = os.path.join(path, name).split('.')[0]
            print(file_dirname)
            command = fr"python ../matlab2python/matlab2python.py {file_dirname}.m -o {file_dirname}.py"
            p = subprocess.run(command)