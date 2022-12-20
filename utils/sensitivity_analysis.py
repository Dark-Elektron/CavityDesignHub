import json

import plotly.express as px
import numpy as np
from icecream import ic
from numpy.random import normal
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import plotly
import pandas as pd
import plotly.graph_objs as go
import itertools
from itertools import combinations
import matplotlib.ticker as mticker
from sympy import *
from sklearn import linear_model


def plot_settings():
    import matplotlib as mpl
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20

    mpl.rcParams['axes.labelsize'] = 24
    mpl.rcParams['axes.titlesize'] = 24
    mpl.rcParams['legend.fontsize'] = 14
    mpl.rcParams['legend.title_fontsize'] = 14

    mpl.rcParams['figure.figsize'] = [10, 6]
    mpl.rcParams['figure.dpi'] = 100


plot_settings()


def P(n, x, symbol=None):
    if symbol:
        sym = symbols(symbol)  # sympy symbol

        if n == 0:
            return 1, sym

        if n == 1:
            return sym, sym

        if n == 2:
            return 0.5 * (3 * sym ** 2 - 1), sym

        if n == 3:
            return 0.5 * (5 * sym ** 3 - 3 * sym), sym

        if n == 4:
            return 0.125 * (35 * sym ** 4 - 30 * sym ** 2 + 3), sym
    else:
        shift = (max(x) + min(x))/2
        scale = (max(x) - min(x))/2

        if n == 0:
            return np.ones(len(x))

        if n == 1:
            x = (np.array(x) - shift) / scale
            return x

        if n == 2:
            x = (np.array(x) - shift) / scale
            return 0.5 * (3 * x ** 2 - 1)

        if n == 3:
            x = (np.array(x) - shift) / scale
            return 0.5 * (5 * x ** 3 - 3 * x)

        if n == 4:
            x = (np.array(x) - shift) / scale
            return 0.125 * (35 * x ** 4 - 30 * x ** 2 + 3)


def He(n, x):
    shift = (max(x) + min(x))/2
    scale = (max(x) - min(x))/2
    if n == 0:
        return np.ones(len(x))

    if n == 1:
        x = (np.array(x) - shift) / scale
        return 2*x

    if n == 2:
        x = (np.array(x) - shift) / scale
        return 4 * x ** 2 - 2

    if n == 3:
        x = (np.array(x) - shift) / scale
        return 8 * x ** 3 - 12 * x

    if n == 4:
        x = (np.array(x) - shift) / scale
        return 16 * x ** 4 - 48 * x ** 2 + 12


def sub_lists(l):
    # initializing empty list
    comb = []

    # Iterating till length of list
    for i in range(len(l) + 1):
        # Generating sub list
        comb += [list(j) for j in combinations(l, i)]
    # Returning list
    return comb[1:-1]


def sobol(d, obj, df_=None):
    if isinstance(obj, list):
        symbols_dict, poly = obj

    ic("Sobol")
    dd = {}
    for k, v in d.items():
        dd[k] = np.linspace(min(v), max(v), 2)

    mg = np.meshgrid(*dd.values())

    df_ = pd.DataFrame()
    x = {}
    for i, key in enumerate(d.keys()):
        vv = mg[i].flatten()
        df_[key] = vv
        x[key] = vv

    for key, ob in obj.items():
        df_[key] = eval(ob)

    # ic(df_.loc[0])
    # vvv = []
    # r_v = ['lh_1', 'lh_3', 'lh_4', 'dh_3', 'alpha_h', 'ch_2', 'r_{cyl}', 'offset_y']
    # for k, v in df_.iterrows():
    #     symbols_dict = update_symbols_dict(symbols_dict, v[r_v])
    #     # print(symbols_dict)
    #     pv = poly.subs(symbols_dict)
    #     vvv.append(pv)
        # print(pv)

    # ic(vvv)

    # get subsets
    power_set = sub_lists(d.keys())
    # ic(power_set)

    Sj = {}
    STi = {}
    for ob in obj.keys():
        Sj[ob] = {}
        STi[ob] = {}
        for subset in power_set:
            # ic(len(power_set))
            # ic(power_set)
            var = np.var(df_.groupby(subset)[ob].mean()) / np.var(df_[ob])

            if len(subset) == 1:
                # print(subset, np.var(df_.groupby(subset)[ob].mean()), np.var(df_[ob]))
                # if subset == ['a2']:
                #     print(df_.groupby(subset)[ob].mean())
                #     print(df_[ob])
                Sj[ob][f'V[E[{ob}|{",".join(map(str, subset))}]]'] = var

                # total sobol
                subset_inv = [s for s in d.keys() if s not in subset]
                var_i = 1 - np.var(df_.groupby(subset_inv)[ob].mean())/np.var(df_[ob])
                STi[ob][f'V[E[{ob}|{",".join(map(str, subset))}]]'] = var_i
            else:
                power_subset = sub_lists(subset)
                for rand_var in power_subset:
                    var = var - Sj[ob][f'V[E[{ob}|{",".join(map(str, rand_var))}]]']

                Sj[ob][f'V[E[{ob}|{",".join(map(str, subset))}]]'] = var

    return Sj, STi


def sobol_df(x, obj, df_):
    power_set = sub_lists(x)

    Sj = {}
    STi = {}
    for ob in obj:
        Sj[ob] = {}
        STi[ob] = {}
        for subset in power_set:
            # ic(len(power_set))
            # ic(power_set)
            # ic(df_[ob], np.var(df_[ob]))
            var = np.var(df_.groupby(subset)[ob].mean()) / np.var(df_[ob])

            if len(subset) == 1:
                # print(subset, np.var(df_.groupby(subset)[ob].mean()), np.var(df_[ob]))
                # if subset == ['a2']:
                #     print(df_.groupby(subset)[ob].mean())
                #     print(df_[ob])
                Sj[ob][f'V[E[{ob}|{",".join(map(str, subset))}]]'] = var

                # total sobol
                subset_inv = [s for s in x if s not in subset]
                var_i = 1 - np.var(df_.groupby(subset_inv)[ob].mean())/np.var(df_[ob])
                STi[ob][f'V[E[{ob}|{",".join(map(str, subset))}]]'] = var_i
            else:
                power_subset = sub_lists(subset)
                for rand_var in power_subset:
                    var = var - Sj[ob][f'V[E[{ob}|{",".join(map(str, rand_var))}]]']

                Sj[ob][f'V[E[{ob}|{",".join(map(str, subset))}]]'] = var

    return Sj, STi


def pce():
    pass


def check_distribution(df):
    v = ['A', 'B', 'a2', 'b3', 'Ri', 'Req', 'Epk/Eacc', 'Bpk/Eacc', 'R/Q']
    # v = ['A', 'B', 'a', 'b', 'Ri', 'Req', 'Epk/Eacc_1', 'Bpk/Eacc_1', 'Rsh/Q_1']
    # v = ['Epk/Eacc', 'Bpk/Eacc', 'R/Q']

    for p in v:
        # df[p] = df[p] / df[p].abs().max()
        # print(df['B'].describe())
        df[p].plot.density()


def update_symbols_dict(sym_dict, values):
    assert len(sym_dict) == len(values)
    new_sym_list = {}
    for i, a in enumerate(sym_dict.keys()):
        new_sym_list[a] = values[i]

    return new_sym_list


def projection():

    return c, p


def regression(df, poly_list_sym, obj, symbols_dict):
    reg = linear_model.LinearRegression(fit_intercept=True)
    # build matrix
    A = []
    # print(poly_list_sym)
    for k, v in df.iterrows():
        symbols_dict = update_symbols_dict(symbols_dict, v[r_v])
        aa = [p.subs(symbols_dict) if not isinstance(p, np.int32) else 1 for p in poly_list_sym.values()]
        A.append(aa)

    poly, coef_dict = {}, {}
    for ob in obj:
        b = df[ob].to_numpy(dtype='float')
        reg.fit(np.array(A, dtype='float'), b.reshape(-1, 1))

        coef_dict[ob] = reg.coef_.copy()[0]
        coef_dict[ob][0] = reg.intercept_[0]
        ic(reg.score(A, b))
        ic(x)

        poly[ob] = np.sum(np.array(list(poly_list_sym.values()))*coef_dict[ob])

        symbols_dict = update_symbols_dict(symbols_dict, df.loc[1, r_v])
        ic(poly[ob].subs(symbols_dict))

    return poly, coef_dict


def generate_sobol_sequence(dim, index, columns, bounds):
    from scipy.stats import qmc
    sampler = qmc.Sobol(d=dim, scramble=False)
    sample = sampler.random_base2(m=index)
    ic(qmc.discrepancy(sample))
    u_bounds, l_bounds = bounds.T

    sample = qmc.scale(sample, l_bounds, u_bounds)
    # print(sample)

    df = pd.DataFrame(sample, columns=columns)
    # ic(df)
    # writePath = r'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\dqw_random_vector.txt'
    #
    # with open(writePath, 'w') as f:
    #     dfAsString = df.to_string(header=True, index=False)
    #     f.write(dfAsString)
    return df


def simplySupportedBeam():

    var = ['b', 'h', 'L', 'E', 'p']
    var_int = np.array([[0.15, 0.0075], [0.3, 0.015], [5, 0.05], [3e10, 4.5e9], [1e4, 2e3]])  # big interval
    n = 1e1

    # df = generate_sobol_sequence(5, 8, var, var_int)

    dd = {}
    for k, v in enumerate(var):
        dd[v] = np.linspace(min(var_int[k]), max(var_int[k]), int(n))
    mg = np.meshgrid(*dd.values())
    df = pd.DataFrame()
    for i, rand_var in enumerate(var):
        df[rand_var] = mg[i].flatten()

    print(df)

    b, h, L, E, p = df.to_numpy().T
    print(b)

    # The beam is considered primatic, therefore:
    I = b * h**3 / 12  # the moment of inertia

    # now for the actual execution we use a vectorized formula:
    Y = np.zeros((len(L), 9))
    for jj in range(1, 10):
        # calculate the xi values:
        xi = jj/10 * L

        Y[:, jj-1] = -p * xi * (L**3 - 2*xi**2 * L + xi**3)/(24*E*I)

    Y = pd.concat([df, pd.DataFrame(Y, columns=['Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8', 'Y9'])], axis=1)
    return Y


if __name__ == '__main__':

    filename = fr'D:\Dropbox\CEMCodesHub\C800MHz\PostprocessingData\Data\GridSimulation_Data.xlsx'
    filename = fr'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\HC_Smax_Fmax_Data_7var_2p.xlsx'
    filename = fr'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\DQW_Smax_Fmax_Data_11var_2p.xlsx'
    df = pd.read_excel(filename, 'Sheet1')
    # ic(df)

    # random_var_order = [5, 5, 5, 5, 5]
    # random_var_order = [7, 7, 7, 7, 7]
    p_order, truncation = 2, 1

    # r_v = ['lh_1', 'lh_3', 'lh_4', 'dh_3', 'alpha_h', 'ch_2', 'r_{cyl}', 'offset_y']
    r_v = ['shaft_y', 'bend_out_sec_prob_length', 'bend_out_cap_gap', 'Window_margin', 'Shift_from_center',
           'Cap_y', 'Cap_thickness', 'Cap_Height', 'Cap_Gap', 'Bend_out_chamber_Length', 'BP_HOM_Penetration']

    random_var_order = [p_order for v in r_v]
    # random_var_order = [1, 1, 1, 1]
    rvo = [[i for i in range(x+1)] for x in random_var_order]
    ic(rvo)
    alpha = list(itertools.product(*rvo))
    x = {}
    for r in r_v:
        x[r] = df[r]
    # ic(x)

    # x = {'A': df.A, 'B': df.B,
    #      'a2': df.a2, 'b3': df.b3,
    #      'Ri': df.Ri
    #      }

    # check_distribution(df)
    # plt.legend()
    # plt.show()

    # z = (1-x[0])**2 + 100*(x[1] - x[0]**2)**2

    poly_list = {}
    poly_chaos_ex = {}
    ic(len(alpha))
    poly_list_sym = {}
    symbols_dict = {}
    for a in alpha:
        # truncate
        if sum(a) <= truncation:
        # if np.linalg.norm(a) <= truncation:
            poly_list[f'{a}'] = [P(i, x[list(x.keys())[j]]) for j, i in enumerate(a)]

            ll = []
            for j, i in enumerate(a):
                pc, sy = P(i, x[list(x.keys())[j]], r_v[j])
                ll.append(pc)

                if len(symbols_dict) <= len(a):
                    symbols_dict[sy] = 0

            poly_list_sym[f'{a}'] = np.prod(ll)

            poly_chaos_ex[f'{a}'] = [f"P({i}, x['{list(x.keys())[j]}'])" for j, i in enumerate(a)]

    print(symbols_dict)
    ic(poly_chaos_ex)
    ic(len(poly_chaos_ex))
    # ic(poly_list_sym)

    # obj = ['Epk/Eacc', 'Bpk/Eacc', 'R/Q']
    obj = ['S_0', 'f_0', 'BW']

    # update symbols_dict with values
    # ic(df)
    for i in range(len(df)):
        symbols_dict = update_symbols_dict(symbols_dict, df.loc[i, r_v])

    poly, coef = regression(df, poly_list_sym, obj, symbols_dict)

    # obj = ['f_0']
    c = {}
    poly_sym_d = {}
    ic(len(poly_list), len(coef))
    for ob in obj:
        c[ob] = {}
        poly_sym = 0
        ci = 0
        for k, v in poly_list.items():
            p = np.prod(np.vstack(v), axis=0)

            c[ob][k] = np.dot(df[ob], p)/np.dot(p, p)
            poly_sym += c[ob][k]*poly_list_sym[k]
            # poly_sym += coef[ob][ci] * poly_list_sym[k]
            ci += 1

        poly_sym_d[ob] = poly_sym

    ic(c['f_0'])
    ic(poly_sym_d['f_0'])
    # ic(np.mean(df['Epk/Eacc']))

    # free some memory
    poly_list = None

    # build polynomials
    pp = {}
    obj_pce = {}
    for ob in obj:
        obj_pce[ob] = ' + '.join(map(str, [f"{coeff}*{'*'.join(map(str, poly_chaos_ex[key]))}" for key, coeff in c[ob].items()]))
        # obj[ob] = pp[ob]
    # rand_var = ['x1', 'x2', 'x3', 'x4']
    # for o
    # obj = {'f1': pp['Epk/Eacc'], 'f2': pp['Bpk/Eacc'], 'f3': pp['R/Q']}

    # check accuracy of polynomial

    # df = simplySupportedBeam()
    # x = ['b', 'h', 'L', 'E', 'p']
    # obj = ['Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8', 'Y9']
    # r_v = x
    # ic(df)
    ic(obj_pce['f_0'])
    # Sj, STi = sobol_df(x, obj, df)
    Sj, STi = sobol(x, obj_pce, df)

    # save sobol indices
    with open(fr"D:\Dropbox\CEMCodesHub\utils\Sobol\Sj_p{p_order}_t{truncation}.json", 'w') as file:
        file.write(json.dumps(Sj, indent=4, separators=(',', ': ')))

    with open(fr"D:\Dropbox\CEMCodesHub\utils\Sobol\STi_p{p_order}_t{truncation}.json", 'w') as file:
        file.write(json.dumps(STi, indent=4, separators=(',', ': ')))

    # ic(Sj)

    for o in obj:
        ic(sum(Sj[o].values()))

    # plot
    fig, ax = plt.subplots(1, 2)
    width = 0.5

    # for ob in obj:
    # rand_var_dict = {0: 'A', 1: 'B', 2: 'a', 3: 'b', 4: 'Ri'}
    # rand_var_dict = {0: 'lh1', 1: 'lh3', 2: 'lh4', 3: 'dh3', 4: 'alpha_h', 5: 'ch2', 6: 'r_cyl', 7: 'offset_y'}
    rand_var_dict = {}
    for i, v in enumerate(x):
        rand_var_dict[i] = v

    bottom, bottom1 = np.zeros(len(obj)), np.zeros(len(obj))
    for key in x:
    # for key in x.keys():
        if key == 0:
            ax[0].bar(obj, [Sj[ob][f'V[E[{ob}|{key}]]'] for ob in obj], width, label=r"$\mathbf{" + key + '}$')
            ax[1].bar(obj, [STi[ob][f'V[E[{ob}|{key}]]'] for ob in obj], width, label=r"$\mathbf{" + key + '}$')
        else:
            ax[0].bar(obj, [Sj[ob][f'V[E[{ob}|{key}]]'] for ob in obj], width, label=r"$\mathbf{" + key + '}$', bottom=bottom)
            ax[1].bar(obj, [STi[ob][f'V[E[{ob}|{key}]]'] for ob in obj], width, label=r"$\mathbf{" + key + '}$', bottom=bottom1)

        bottom += np.array([Sj[ob][f'V[E[{ob}|{key}]]'] for ob in obj])
        bottom1 += np.array([STi[ob][f'V[E[{ob}|{key}]]'] for ob in obj])

    # ylabel = ['$S_j$', '$ST_i$']
    ylabel = ['$S_j$', '$ST_i$']
    for i, a in enumerate(ax):
        ticks_loc = a.get_xticks()
        a.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        # a.set_xticklabels(['Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8', 'Y9'])
        a.set_xticklabels(obj)
        a.set_ylabel(ylabel[i])

    # Add a table at the bottom of the axes
    S = [Sj, STi]
    for i, s in enumerate(S):
        cellText = []
        rowLabels = []
        colLabels = []
        n = 0
        for k, v in s.items():
            rowLabels.append(r"$\mathbf{" + k + '}$')

            temp_col, temp_cellText = [], []
            for xx, xv in zip(list(v.keys()), list(v.values())):
                if xx in [f'V[E[{k}|{p}]]' for p in r_v]:
                    temp_col.append(r"$\mathbf{" + f"{xx.replace(k, '.')}" + '}$')
                    temp_cellText.append(r"$\mathbf{" + f"{round(xv, 4)}" + '}$')

            if len(temp_col) != 0:
                if n == 0:
                    colLabels = temp_col
                cellText.append(temp_cellText)

            n += 1

        # vals = cellText
        # print(cellText)
        # # norm = plt.Normalize(np.min(vals) - 1, np.max(vals) + 1)
        # # print(norm)
        # colours = plt.cm.binary(normal(vals), 0.5)
        rcolors = plt.cm.BuPu(np.full(len(rowLabels), 0.5))
        ccolors = plt.cm.BuPu(np.full(len(temp_col), 0.5))
        # print(rowLabels, colLabels)
        the_table = ax[i].table(cellText=cellText,
                                   rowLabels=rowLabels,
                                   colLabels=colLabels,
                                   rowColours=rcolors,
                                   colColours=ccolors,
                                   # cellColours=colours,
                                   loc='bottom', bbox=[0.0, -0.35, 1, .28])
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(8)

    plt.subplots_adjust(bottom=0.3)
    plt.subplots_adjust(top=0.93)

    lines, labels = ax[0].get_legend_handles_labels()

    fig.legend(lines, labels, loc="upper center", ncol=int(len(rand_var_dict.keys())), fancybox=True, shadow=False, bbox_to_anchor=(0.5, 1.0))

    # ax[0].legend(loc="upper center", ncol=int(len(rand_var_dict.keys())), bbox_to_anchor=(0.0, 1.1),
    #        fancybox=True, shadow=True)
    # ax[1].legend(loc="upper center", ncol=int(len(rand_var_dict.keys())), bbox_to_anchor=(0.0, 1.1),
    #           fancybox=True, shadow=True)

    plt.tight_layout()
    plt.show()
