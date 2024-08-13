import os.path
import sys
import pandas as pd
import numpy as np


def read_output_from_cst_sweep(sim_folder, folders, requests):
    df_ = pd.DataFrame()
    for request in requests:
        d = pd.DataFrame()
        for folder in folders:
            d = pd.concat([d, pd.read_csv(fr"{sim_folder}\{folder}\Export\{request}.txt", sep="\t", header=None)])
            # print(d)

        df_[request] = d.loc[:, 1]

    return df_


if len(sys.argv) != 6:
    print("Usage: python myscript.py param_file output_file partitions")
    sys.exit(1)

param_file = sys.argv[1]
output_file = sys.argv[2]
partitions = sys.argv[3]
nodes_only = sys.argv[4]
num_responses = sys.argv[5]

assert isinstance(partitions, int) or isinstance(partitions, float), print("Partitions argument must be integer.")
partitions = int(partitions)

# print(f"parameter_file: {param_file}, output_file: {output_file}")

# open parameter file and get parameters
df = pd.read_csv(param_file, sep='\s+', header=None)
num_in_vars = int(df.loc[0, 0])

pars_in = df.loc[1:num_in_vars, 0]

if nodes_only != 'ONLY_NODES':
    # process output from cst run
    current_dir = os.getcwd()
    # folders = [f"{os.getcwd()}/model_{partition+1}" for partition in range(int(partitions))]
    # req = ["Frequency (Multiple Modes)", "Q-Factor (lossy E) (Multiple Modes)", "R over Q beta=1 (Multiple Modes)",
    #        "RQT", "RQT_kOhm_m", "Z_kOhm", "Z_T_kOhm_m"]
    #
    # df = read_output_from_cst_sweep(current_dir, folders, requests)
    filename = fr'{current_dir}\cubature_nodes.xlsx'
    df_data = pd.read_excel(fr"{filename}", 'Sheet1', engine='openpyxl')

    for indx, row in df_data.iterrows():
        tolerance = 1e-3
        if np.allclose(np.around(row.tolist()[:num_in_vars], 3), np.around(pars_in.tolist(), 3), rtol=tolerance, atol=tolerance):
            # write output
            out = row.tolist()[num_in_vars:]
            with open(output_file, 'w') as f:
                for o in out:
                    f.write(f"{o:20.10e}     f\n")
else:
    out = np.ones(int(num_responses))
    with open(output_file, 'w') as f:
        for o in out:
            f.write(f"{o:20.10e}     f\n")
