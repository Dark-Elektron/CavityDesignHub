import random

import numpy as np
import oapackage
import pandas
import pandas as pd
from icecream import ic


def pareto_front(x, y):
    def reverse_list(l, goal):
        if goal == 'max':
            return l  # to find the pareto maxima
        else:
            return [-x for x in l]  # to find the pareto minima

    datapoints = np.array([reverse_list(x, 'min'), reverse_list(y, 'min')])

    pareto = oapackage.ParetoDoubleLong()

    for ii in range(0, datapoints.shape[1]):
        w = oapackage.doubleVector((datapoints[0, ii], datapoints[1, ii]))
        pareto.addvalue(w, ii)
    # pareto.show(verbose=1)  # Prints out the results from pareto

    lst = pareto.allindices()  # the indices of the Pareto optimal designs

    optimal_datapoints = datapoints[:, lst]

    return reverse_list(optimal_datapoints[0, :], 'min'), reverse_list(optimal_datapoints[1, :], 'min')


def crossover(df, generation, f):  # , rq, grq
    df_co = pd.DataFrame(columns=['key', 'x', 'y'])
    for i in range(f):
        # (<obj>[<rank>][<variable>] -> (b[c[1]][0]

        df_co.loc[i] = [f"G{generation}_C{i}_CO",
                     df.loc[np.random.randint(5)]["x"],  # A
                     df.loc[np.random.randint(5)]["y"]]  # B]
    return df_co


def mutation(df, n):
    df_ng_mut = pd.DataFrame(columns=['key', 'x', 'y'])
    df_ng_mut.loc[:, 'A'] = df.loc[:, "A"] * random.uniform(0.85, 1.5)
    df_ng_mut.loc[:, 'B'] = df.loc[:, "B"] * random.uniform(0.85, 1.5)

    key1, key2 = [], []
    for i in range(len(df_ng_mut)):
        key1.append(f"G{n}_C{i}_M")

    df_ng_mut.loc[:, 'key'] = key1

    return df_ng_mut


def chaos(f, n):
        data = {'key': [f"G{n}_C{i}_CH" for i in range(f)],
                'x': random.sample(list(np.linspace(0.5, 1.5, f * 10)), f),
                'y': random.sample(list(np.linspace(0.5, 1.5, f * 10)), f)}

        df = pd.DataFrame.from_dict(data)
        return df


def ea(self, df, obj, constraints, n):
    # optimize by page rank
    # remove the lowest ranking members

    # compare with global dict and remove duplicates
    compared_cols = ['x', 'y']
    if not self.df_global.empty:
        df = df.loc[~df.set_index(compared_cols).index.isin(df_global.set_index(compared_cols).index)]  # this line of code removes duplicates

    # replace Req with tuned Req
    df[['f1', 'f2']] = [(x-1)**2 + 3*(y-1)**2, 4*x**2 + y**2 + x*y]

    ######################################################
    # filter by constraints
    for const in self.constraints:
        c = const.split(" ")

        if c[1] == '>':
            df = df.loc[(df[f'{c[0]}'] > float(c[2]))]
        elif c[1] == '<':
            df = df.loc[(df[f'{c[0]}'] < float(c[2]))]
        elif c[1] == '<=':
            df = df.loc[(df[f'{c[0]}'] <= float(c[2]))]
        elif c[1] == '>=':
            df = df.loc[(df[f'{c[0]}'] >= float(c[2]))]
        elif c[1] == '==':
            df = df.loc[(df[f'{c[0]}'] == float(c[2]))]

    # update with global dataframe
    if not self.df_global.empty:
        df = pd.concat([self.df_global, df], ignore_index=True)

    # reset total rank
    df['total_rank'] = 0

    # rank shapes by objectives
    for i, obj in enumerate(self.objectives):

        if obj[0] == "min":
            df[f'rank_{obj[1]}'] = df[obj[1]].rank() * self.weights[i]
        elif obj[0] == "max":
            df[f'rank_{obj[1]}'] = df[obj[1]].rank(ascending=False) * self.weights[i]
        elif obj[0] == "equal":  # define properly later
            continue

        # if 'total_rank' in df.columns:
        df[f'total_rank'] = df[f'total_rank'] + df[f'rank_{obj[1]}']

    # reorder
    ic(df)
    tot = df.pop(f'total_rank')
    ic(tot)
    df[f'total_rank'] = tot / sum(self.weights)  # normalize by sum of weights
    ic(df)

    # order shapes by rank
    df = df.sort_values(by=['total_rank'])
    df = df.reset_index(drop=True)

    # pareto condition
    reorder_indx = self.pareto_front(df)
    df = df.loc[reorder_indx, :]
    df = df.reset_index(drop=True)
    ic(df)

    # update global
    if len(df) > self.tuneUI.sb_Max_Table_Size.value():
        # self.df_global = df.loc[0:self.tuneUI.cb_Max_Table_Size.value(), :]
        df_global = df
    else:
        df_global = df
    ic(df_global)

    # save dataframe
    filename = fr"{self.projectDir}\SimulationData\SLANS\Generation{n}.xlsx"
    recursive_save(self.df_global, filename, reorder_indx)

    # crossover
    print("Crossover")
    df_cross = crossover(df, n, self.tuneUI.sb_Crossover_Factor.value())  # , elites["GR/Q
    # ic(df_cross)

    # mutation
    print("Mutation")
    df_mutation = mutation(df, n)
    # ic(df_mutation)

    # chaos
    print("Chaos")
    df_chaos = chaos(self.tuneUI.sb_Chaos_Factor.value(), n)
    # ic(df_chaos)

    # take elites from previous generation over to next generation
    df_ng = pd.concat([df_cross, df_mutation, df_chaos], ignore_index=True)

    # update dictionary
    df = df_ng

    n += 1
    print(n)
    print("=" * 80)
    if n < self.ng_max:
        return self.ea(n)
    else:
        return


f1_arr, f2_arr = [], []
for x in np.linspace(-0.5, 1.5, 1000):
    for y in np.linspace(-0.5, 1.5, 100):
        f1 = (x-1)**2 + 3*(y-1)**2
        f2 = 4*x**2 + y**2 + x*y

        f1_arr.append(f1)
        f2_arr.append(f2)


f1_p, f2_p = pareto_front(f1_arr, f2_arr)

import matplotlib.pyplot as plt
plt.scatter(f1_p, f2_p, s=5)
plt.show()
