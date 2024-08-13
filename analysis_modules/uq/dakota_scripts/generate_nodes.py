import subprocess
import os
import sys

import numpy as np
import pandas as pd


# Get current path
sCurPath = os.path.abspath(".")

if len(sys.argv) != 3:
    print("Usage: python myscript.py param_file output_file partitions")
    sys.exit(1)

param_file = sys.argv[1]
output_file = sys.argv[2]

# open parameter file and get parameters
df = pd.read_csv(param_file, sep='\s+', header=None)
num_in_vars = int(df.loc[0, 0])

pars_in = df.loc[1:num_in_vars, 0]

out = np.ones(1)
with open(output_file, 'w') as f:
    for o in out:
        f.write(f"{o:20.10e}     f\n")

