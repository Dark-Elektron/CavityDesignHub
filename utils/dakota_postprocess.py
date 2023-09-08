import csv
import re

import pandas as pd
from icecream import ic
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker


def plot_sobol_indices(filepath, objectives, which=None, kind='stacked', orientation='vertical', normalise=True):
    if which is None:
        which = ['Main']

    start_keyword = "Main"
    interaction_start_keyword = "Interaction"
    pattern = r'\s+(-?\d\.\d+e[+-]\d+)\s+(-?\d\.\d+e[+-]\d+)\s+(\w+)'
    pattern_interaction = r'\s*(-?\d+\.\d+e[-+]\d+)\s+(\w+)\s+(\w+)\s*'

    with open(filepath, "r") as file:
        # read the file line by line
        lines = file.readlines()

        # initialize a flag to indicate when to start and stop recording lines
        record = False
        record_interaction = False

        # initialize a list to store the lines between the keywords
        result = {}
        result_interaction = {}
        count = 0
        # loop through each line in the file
        for line in lines:
            # check if the line contains the start keyword
            if start_keyword in line:
                # if it does, set the flag to start recording lines
                record = True
                result[count] = []
                continue

            if interaction_start_keyword in line:
                record_interaction = True
                result_interaction[count] = []
                continue

            # if the flag is set to record, add the line to the result list
            if record:
                if re.match(pattern, line):
                    result[count].append(re.findall("\S+", line))
                else:
                    record = False
                    count += 1

            if record_interaction:
                if re.match(pattern_interaction, line):
                    result_interaction[count-1].append(re.findall("\S+", line))
                else:
                    record_interaction = False

    ic(result, type(result))
    ic(result_interaction)
    # result = result_interaction
    df_merge = pd.DataFrame(columns=['main', 'total', 'vars'])
    ic(df_merge)
    # print the lines between the keywords
    # df_merge_list = []
    for i, (k, v) in enumerate(result.items()):
        df = pd.DataFrame(v, columns=['main', 'total', 'vars'])
        df = df.astype({'main': 'float', 'total': 'float'})
        # df_merge_list.append(df)
        # ic(df)
        # ic(df_merge)
        if i == 0:
            df_merge = df
        else:
            df_merge = pd.merge(df_merge, df, on='vars', suffixes=(f'{i}', f'{i+1}'))
        # df.plot.bar(x='var', y='Main')
    # ic(df_merge_list)
    # df_merge = pd.merge(df_merge_list, on='vars')
    # ic(df_merge)

    df_merge_interaction = pd.DataFrame(columns=['interaction', 'var1', 'var2'])
    for i, (k, v) in enumerate(result_interaction.items()):
        df = pd.DataFrame(v, columns=['interaction', 'var1', 'var2'])
        df = df.astype({'interaction': 'float'})
        if i == 0:
            df_merge_interaction = df
        else:
            df_merge_interaction = pd.merge(df_merge_interaction, df, on=['var1', 'var2'])

        # combine var columns
        # df_merge_interaction["vars"] = df_merge_interaction[["var1", "var2"]].agg(','.join, axis=1)
        # df.plot.bar(x='var', y='Main')

    ic(df_merge)
    ic(df_merge_interaction)

    if normalise:
        # normalise dataframe columns
        for column in df_merge.columns:
            if 'main' in column or 'total' in column:
                df_merge[column] = df_merge[column].abs() / df_merge[column].abs().sum()

    ic(df_merge)
    # filter df
    for w in which:
        if w.lower() == 'main' or w.lower() == 'total':
            dff = df_merge.filter(regex=f'{w.lower()}|vars')
        else:
            # create new column which is a combination of the two variable names
            dff = df_merge_interaction.filter(regex=f'{w.lower()}|vars')
            if not dff.empty:
                dff['vars'] = df_merge_interaction[['var1', 'var2']].apply(lambda x: '_'.join(x), axis=1)
        ic(dff)

        cmap = 'tab20'

        if not dff.empty:
            if kind.lower() == 'stacked':
                dff_T = dff.set_index('vars').T
                if orientation == 'vertical':
                    ax = dff_T.plot.bar(stacked=True, rot=0, cmap=cmap)
                    ax.set_xlim(left=0)
                    plt.legend(bbox_to_anchor=(1.04, 1), ncol=2)
                else:
                    ax = dff_T.plot.barh(stacked=True, rot=0, cmap=cmap, edgecolor='k')
                    ax.set_xlim(left=0)
                    ax.set_yticklabels(objectives)
                    plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=7, loc='lower left', mode='expand')
            else:
                if orientation == 'vertical':
                    ax = dff.plot.bar(x='vars', stacked=True, cmap=cmap)
                    ax.set_xlim(left=0)
                    ax.axhline(0.03, c='k')
                    plt.legend(bbox_to_anchor=(1.04, 1), ncol=2)
                else:
                    ax = dff.plot.barh(x='vars', stacked=True, cmap=cmap)
                    ax.set_xlim(left=0)
                    ax.axvline(0.03, c='k')
                    plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=7, loc='lower left', mode='expand')

            plt.tight_layout()
            plt.show()
        else:
            ic(f"No {w} found.")


def quadrature_nodes_to_cst_par_input(filefolder, n=2):
    filepath = fr"{filefolder}\sim_result_table.dat"
    df = pd.read_csv(filepath, sep='\s+')

    # delete unnecessary columns
    df.drop(df.filter(regex='response|interface|eval_id').columns, axis=1, inplace=True)
    df.to_excel(fr"{filefolder}\cubature_nodes_pars.xlsx", index=False)

    # save parts
    row_partition = len(df.index)//n
    for i in range(n):
        if i < n-1:
            df_part = df.loc[i*row_partition:(i+1)*row_partition-1]
        else:
            df_part = df.loc[i * row_partition:]

        df_part.to_csv(fr"{filefolder}\cst_par_in_{i+1}.txt", sep="\t", index=None)


def get_pce(filefolder):

    filepath = fr"{filefolder}\uq_pce_expansion.dat"
    df = pd.read_csv(filepath, sep='\s+', header=None)
    ic(df)
    poly = 0
    for row in df.iterrows():
        poly += 0


def quote_to_float(value):
    if isinstance(value, str) and value.startswith('"') and value.endswith('"'):
        return float(value[1:-1])
    else:
        return value


def combine_params_output(folder, N):
    for i in range(N):
        if i == 0:
            df = pd.read_csv(f'{folder}/m{(i+1):02d}.csv', engine='python', skipfooter=1)
        else:
            df = pd.concat([df, pd.read_csv(f'{folder}/m{(i+1):02d}.csv', engine='python', skipfooter=1)])

    # rearrange column according to the reference column order
    df_reference = pd.read_excel(fr"{folder}\cubature_nodes_pars.xlsx")
    columns = list(df_reference.columns)
    # check if 3D Run ID in column and drop if yes
    if ' 3D Run ID' in list(df.columns):
        columns.append(' 3D Run ID')

    columns = list(df_reference.columns) + (df.columns.drop(columns).tolist())
    ic(columns)

    df = df[columns]
    ic(df)

    df.to_excel(fr"{folder}\cubature_nodes.xlsx", index=False)


if __name__ == '__main__':
    plt.rcParams["figure.figsize"] = (6, 2)
    filefolder = fr"C:\Users\sosoho\DakotaProjects\COMPUMAG\ConferenceResults\HC_MC_10%"
    #
    # obj = [r"$Q_\mathrm{ext, FM}$", r"$\max(Q_\mathrm{ext, dip})$"]
    # obj = [r"$f(S_\mathrm{max})~\mathrm{[MHz]}$", r"$f(S_\mathrm{min})~\mathrm{[MHz]}$"]
    obj = [r"$S_\mathrm{max}~\mathrm{[dB]}$", r"$S_\mathrm{min}~\mathrm{[dB]}$", r"$f(S_\mathrm{max})~\mathrm{[MHz]}$", r"$f(S_\mathrm{min})~\mathrm{[MHz]}$"]
    plot_sobol_indices(fr"{filefolder}\dakota_HC.out", obj, ['main', 'Total', 'Interaction'], kind='stacked', orientation='horizontal')#, normalise=False

    # quadrature_nodes_to_cst_par_input(filefolder, n=10)
    # combine_params_output(filefolder, 10)

    # get_pce(filefolder)

