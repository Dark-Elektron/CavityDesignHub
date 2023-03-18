import re

import pandas as pd
from icecream import ic
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker


def plot_sobol_indices(filepath, objectives, which='Main', kind='stacked', orientation='vertical'):
    start_keyword = "Main"
    end_keyword = "Interaction"

    with open(filepath, "r") as file:
        # read the file line by line
        lines = file.readlines()

        # initialize a flag to indicate when to start and stop recording lines
        record = False

        # initialize a list to store the lines between the keywords
        result = {}
        count = 0
        # loop through each line in the file
        for line in lines:
            # check if the line contains the start keyword
            if start_keyword in line:
                # if it does, set the flag to start recording lines
                record = True
                result[count] = []
                continue

            # check if the line contains the end keyword
            if end_keyword in line:
                # if it does, set the flag to stop recording lines
                record = False
                count += 1

            # if the flag is set to record, add the line to the result list
            if record:
                result[count].append(re.findall("\S+", line))

    df_merge = pd.DataFrame(columns=['Main', 'Total', 'var'])
    # print the lines between the keywords
    for i, (k, v) in enumerate(result.items()):
        df = pd.DataFrame(v, columns=['Main', 'Total', 'var'])
        df = df.astype({'Main': 'float', 'Total': 'float'})
        if i == 0:
            df_merge = df
        else:
            df_merge = pd.merge(df_merge, df, on='var')
        # df.plot.bar(x='var', y='Main')

    ic(df_merge)
    # filter df
    dff = df_merge.filter(regex=f'{which}|var')

    cmap = 'tab20'

    if kind.lower() == 'stacked':
        # filter df
        dff = df_merge.filter(regex=f'{which}|var')
        dff_T = dff.set_index('var').T
        if orientation == 'vertical':
            ax = dff_T.plot.bar(stacked=True, rot=0, cmap=cmap)
            plt.legend(bbox_to_anchor=(1.04, 1), ncol=2)
        else:
            ax = dff_T.plot.barh(stacked=True, rot=0, cmap=cmap, edgecolor='k')
            ax.set_yticklabels(objectives)
            plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=7, loc='lower left', mode='expand')
    else:
        if orientation == 'vertical':
            ax = dff.plot.bar(x='var', stacked=True, cmap=cmap)
            ax.axhline(0.03, c='k')
            plt.legend(bbox_to_anchor=(1.04, 1), ncol=2)
        else:
            ax = dff.plot.barh(x='var', stacked=True, cmap=cmap)
            ax.axvline(0.03, c='k')
            plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=7, loc='lower left', mode='expand')

    plt.tight_layout()
    plt.show()


def quadrature_nodes_to_cst_par_input(filefolder, n=2):
    filepath = fr"{filefolder}\sim_result_table.dat"
    df = pd.read_csv(filepath, sep='\s+')

    # delete unnecessary columns
    df.drop(df.filter(regex='response|interface|eval_id').columns, axis=1, inplace=True)

    # save parts
    row_partition = len(df.index)//n
    for i in range(n):
        if i < n-1:
            df_part = df.loc[i*row_partition:(i+1)*row_partition-1]
        else:
            df_part = df.loc[i * row_partition:]

        df_part.to_csv(fr"{filefolder}\cst_par_in_{i}.txt", sep="\t", index=None)


def get_pce(filefolder):

    filepath = fr"{filefolder}\uq_pce_expansion.dat"
    df = pd.read_csv(filepath, sep='\s+', header=None)
    ic(df)
    poly = 0
    for row in df.iterrows():
        poly += 0


def combine_params_output(params_file, output_file):
    pass


if __name__ == '__main__':
    plt.rcParams["figure.figsize"] = (6, 2.5)
    filefolder = fr"C:\Users\sosoho\DakotaProjects\COMPUMAG\HC_Stroud_3"
    # # filepath = fr"C:\Users\sosoho\DakotaProjects\HookCoupler_COMPUMAG\dakota_HC.out"
    # filepath = fr"C:\Users\sosoho\DakotaProjects\HookCoupler_test\dakota_HC.out"
    # # filepath = fr"C:\Users\sosoho\DakotaProjects\HookCoupler_Shahnam\dakota_HC.out"
    # # filepath = fr"C:\Users\sosoho\DakotaProjects\HookCoupler_COMPUMAG_red\dakota_HC.out"
    # # filepath = fr"C:\Users\sosoho\DakotaProjects\Validate results\HC_Stroud_3\dakota_HC.out"
    # filepath = fr"C:\Users\sosoho\DakotaProjects\Validate results\HC_Stroud_5\dakota_HC.out"
    #
    obj = [r"$Q_\mathrm{ext, FM}$", r"$\max(Q_\mathrm{ext, dip})$"]
    plot_sobol_indices(fr"{filefolder}\dakota_HC.out", obj, 'Main', kind='stacked', orientation='horizontal')
    plot_sobol_indices(fr"{filefolder}\dakota_HC.out", obj, 'Total', kind='stacked', orientation='horizontal')

    # quadrature_nodes_to_cst_par_input(filefolder, n=6)
    # get_pce(filefolder)

