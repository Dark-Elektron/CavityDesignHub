import subprocess
import os
import sys

# Get current path
sCurPath = os.path.abspath(".")

# Get the command line arguments passed by DAKOTA
paramsfile = sys.argv[1]
resultsfile = sys.argv[2]

# Run the Python script and capture the output
cmd = ["python3", os.path.join(sCurPath, "py_dakota.py"), paramsfile, resultsfile, "2", "ALL"]
# cmd = ["python3", os.path.join(sCurPath, "py_dakota.py"), paramsfile, resultsfile, "2", "ONLY_NODES"]
process = subprocess.run(cmd)
# output, _ = process.communicate()

# Print the output in the command prompt
# print("Python Output:", output.decode())
