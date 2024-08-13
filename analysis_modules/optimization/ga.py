import os
import random
import shutil
from distutils import dir_util
import matplotlib as mpl
import numpy as np
import oapackage
import pandas as pd
from icecream import ic
from matplotlib import pyplot as plt
from scipy.interpolate import griddata
from scipy.stats import qmc


def plot_settings():
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12

    mpl.rcParams['axes.labelsize'] = 14
    mpl.rcParams['axes.titlesize'] = 14
    mpl.rcParams['legend.fontsize'] = 10
    mpl.rcParams['legend.title_fontsize'] = 12

    mpl.rcParams['figure.dpi'] = 100

    # Set the desired colormap
    plt.rcParams['axes.prop_cycle'] = plt.cycler('color', plt.cm.Set2.colors)


class GA:
    plt.rcParams["figure.figsize"] = (6, 6)

    def __init__(self):
        self.n_interp = 1000
        self.f2_interp = np.zeros(self.n_interp)
        self.interp_error = []
        self.interp_error_avg = []
        self.ng_max = 50
        self.df_global = None

        self.fig, axs = plt.subplot_mosaic([[1, 3]])

        self.ax1 = axs[1]
        self.ax3 = axs[3]

        # self.x_bounds, self.y_bounds = [-0.5, 1.5], [-0.5, 1.5]
        # self.x_bounds, self.y_bounds = [0.1, 1], [0, 5]
        self.x_bounds, self.y_bounds = [0, 1], [0, 1]
        self.fm, self.fc, self.fch = 50, 50, 50

    def run(self, n, df_ng=None):
        if n == 1 and df_ng is None:
            self.df_global = None
            # generate first generation
            df = self.first_men(300, self.x_bounds, self.y_bounds)

            # evaluate objective
            ob1 = self.f1(df['x'], df['y'])
            ob2 = self.f2(df['x'], df['y'])
            df[['f1', 'f2']] = np.array([ob1.tolist(), ob2.tolist()]).T

            # apply constraint
            df = df.loc[(df['f2'] / (0.858 * np.exp(-0.541 * df['f1'])) >= 1) &
                        (df['f2'] / (0.728 * np.exp(-0.295 * df['f1'])) >= 1)]
        else:
            compared_cols = ['x', 'y']
            if not self.df_global.empty:
                df_ng = df_ng.loc[~df_ng.set_index(compared_cols).index.isin(
                    self.df_global.set_index(compared_cols).index)].copy(
                    deep=True)  # this line of code removes duplicates

            # evaluate objective
            ob1 = self.f1(df_ng['x'], df_ng['y'])
            ob2 = self.f2(df_ng['x'], df_ng['y'])
            df_ng.loc[:, ['f1', 'f2']] = np.array([ob1.tolist(), ob2.tolist()]).T

            # apply constraint
            df_ng = df_ng.loc[(df_ng['f2'] / (0.858 * np.exp(-0.541 * df_ng['f1'])) >= 1) & (
                    df_ng['f2'] / (0.728 * np.exp(-0.295 * df_ng['f1'])) >= 1)]

            df = pd.concat([self.df_global, df_ng], ignore_index=True)

        reorder_indx, poc, pareto_list = self.pareto_front(df)
        pareto_shapes = df.loc[pareto_list, ['f1', 'f2']]
        pareto_var_comb = df.loc[pareto_list, ['x', 'y']]

        # fit polynomial # take care for "sharp" functions, polynomial approximation may not be the best
        # sort pareto shapes
        pareto_shapes_sorted = pareto_shapes.sort_values('f1')
        f1 = np.linspace(min(pareto_shapes['f1']), max(pareto_shapes['f1']), self.n_interp)
        # self.ax1.plot(f1, pareto_poly_fit(f1), '-', c='k')
        f2_interp = np.interp(f1, pareto_shapes_sorted['f1'], pareto_shapes_sorted['f2'])
        # self.ax1.plot(f1, f2_interp, c='k')

        error = np.linalg.norm(f2_interp - self.f2_interp)
        self.f2_interp = f2_interp
        self.interp_error.append(error)
        # print(error, np.average(self.interp_error))
        self.interp_error_avg.append(np.linalg.norm(self.interp_error))

        ###############################################################
        # get exact solution
        # x = pareto_var_comb['x']
        # y = (1 / 6) * (np.sqrt(496 * x ** 2 + 232 * x + 1) - 22 * x + 1)
        # self.ax3.scatter(x, y)
        # pareto_var_comb['y_analytic'] = y
        # MSE = np.sum((pareto_var_comb['y_analytic'] - pareto_var_comb['y']) ** 2) / pareto_var_comb.shape[0]
        # print("MSE: ", MSE)
        # self.ax3.scatter(n, MSE, c='k')

        ################################################################
        # calculate mean square error for averaged by grouping
        pareto_grp_avg = pareto_var_comb.groupby(
            pd.cut(df["x"], np.arange(self.x_bounds[0], self.x_bounds[1] + 0.01, 0.01))).mean()
        pareto_grp_avg_y = pareto_var_comb.groupby(
            pd.cut(df["y"], np.arange(self.y_bounds[0], self.y_bounds[1] + 0.01, 0.01))).mean()
        # # get exact solution
        # x = pareto_grp_avg['x']
        # y = (1 / 6) * (np.sqrt(496 * x ** 2 + 232 * x + 1) - 22 * x + 1)
        #
        # # self.ax2.scatter(x, y)
        #
        # pareto_var_comb['y_analytic'] = y
        # MSE_g = np.sum((pareto_grp_avg['y_analytic'] - pareto_grp_avg['y']) ** 2) / pareto_grp_avg.shape[0]
        # self.ax3.scatter(n, MSE_g, facecolors="None", edgecolors='k', lw=2)

        if n + 1 == self.ng_max:
            self.ax1.scatter(pareto_shapes['f1'], pareto_shapes['f2'], marker='o', s=5,
                             label=f'Pareto Front (Gen. {n - 1})', zorder=2, edgecolor='r')
            self.ax1.set_xlabel('f1')
            self.ax1.set_ylabel('f2')
            # self.ax2.scatter(pareto_var_comb['x'], pareto_var_comb['y'], marker='o', s=7, label=f'Pareto {n - 1}')
            # self.ax2.set_xlabel('x')
            # self.ax2.set_ylabel('y')

            # plot analytic solution
            xrange = np.linspace(*self.x_bounds, 50)
            yrange = np.linspace(*self.y_bounds, 50)
            X, Y = np.meshgrid(xrange, yrange)
            # idx = (Y + 9*X >= 6) & (-Y + 9*X >= 1)
            # X[~idx] = 0
            # Y[~idx] = 0

            # F is one side of the equation, G is the other
            F = (8 * X + Y) * (6 * Y - 6)
            G = (2 * Y + X) * (2 * X - 2)

            # self.ax2.contour(X, Y, (F - G), [0], linewidths=[2], colors=['k'])
            # self.ax2.set_xlim(self.x_bounds)
            # self.ax2.set_ylim(self.y_bounds)

            # plot analytic function
            f1 = self.f1(X, Y)
            f2 = self.f2(X, Y)

            df_analytical = pd.DataFrame({'f1': f1.flatten(), 'f2': f2.flatten()})
            # apply constraint
            df_analytical = df_analytical.loc[
                (df_analytical['f2'] / (0.858 * np.exp(-0.541 * df_analytical['f1'])) >= 1) & (
                        df_analytical['f2'] / (0.728 * np.exp(-0.295 * df_analytical['f1'])) >= 1)]

            self.ax1.scatter(df_analytical['f1'], df_analytical['f2'], s=5, color='r', edgecolor='k', facecolor='none')

            # grouped data average
            # self.ax2.scatter(pareto_grp_avg['x'], pareto_grp_avg['y'], marker='P', s=7, label=f'Pareto {n - 1}')
            # self.ax2.scatter(pareto_grp_avg_y['x'], pareto_grp_avg_y['y'], marker='P', s=7, label=f'Pareto {n - 1}')
            # self.ax3.set_yscale('log')

        # print(poc)

        df = df.loc[reorder_indx, :]
        self.df_global = df
        self.recursive_save(df, "../../utils/ga_res.xlsx", reorder_indx, poc)

        # crossover
        # print("Crossover")
        df_cross = self.crossover(df, n, self.fc)  # , elites["GR/Q
        # ic(df_cross)

        # mutation
        # print("Mutation")
        df_mutation = self.mutation(df, n, self.fm)
        # ic(df_mutation)

        # chaos
        # print("Chaos")
        df_chaos = self.chaos(self.fch, n, self.x_bounds, self.y_bounds)
        # ic(df_chaos)

        # take elites from previous generation over to next generation
        df_ng = pd.concat([df_cross, df_mutation, df_chaos], ignore_index=True)

        # apply constraints
        df_ng = df_ng.loc[(df_ng['x'] <= self.x_bounds[1]) & (df_ng['x'] >= self.x_bounds[0])]
        df_ng = df_ng.loc[(df_ng['y'] <= self.y_bounds[1]) & (df_ng['y'] >= self.y_bounds[0])]
        # df_ng = df_ng.loc[(df_ng['y'] + 9*df_ng['x'] >= 6) & (-df_ng['y'] + 9*df_ng['x'] >= 1)]

        n += 1
        ic(n)
        ic("=" * 80)
        if n < self.ng_max:
            return self.run(n, df_ng)
        else:
            ic(df)
            self.ax3.plot(self.interp_error, marker='P', label=r'$\epsilon$')
            # self.ax3.plot(self.interp_error_avg, marker='X', label=r'$||\mathbf{\epsilon}||_2$')
            self.ax3.set_yscale('log')
            self.ax1.legend()
            self.ax3.legend()
            self.fig.suptitle(r"$\min \left\{f_1 = 4x^2 + y^2 + xy, f_2 = (x - 1)^2 + 3(y - 1)^2 \right\}$")
            return

    def run_sensitivity(self, n, df_ng=None):
        if n == 1 and df_ng is None:
            self.df_global = None
            # generate first generation
            df = self.first_men(300, self.x_bounds, self.y_bounds)

            # evaluate objective
            ob1 = self.f1(df['x'], df['y'])
            ob2 = self.f2(df['x'], df['y'])
            df[['f1', 'f2']] = np.array([ob1.tolist(), ob2.tolist()]).T
        else:
            compared_cols = ['x', 'y']
            if not self.df_global.empty:
                df_ng = df_ng.loc[~df_ng.set_index(compared_cols).index.isin(
                    self.df_global.set_index(compared_cols).index)]  # this line of code removes duplicates

            # evaluate objective
            ob1 = self.f1(df_ng['x'], df_ng['y'])
            ob2 = self.f2(df_ng['x'], df_ng['y'])
            df_ng.loc[:, ['f1', 'f2']] = np.array([ob1.tolist(), ob2.tolist()]).T

            df = pd.concat([self.df_global, df_ng], ignore_index=True)

        reorder_indx, poc, pareto_list = self.pareto_front(df)
        pareto_shapes = df.loc[pareto_list, ['f1', 'f2']]
        pareto_var_comb = df.loc[pareto_list, ['x', 'y']]

        # get exact solution
        x = pareto_var_comb['x']
        y = (1 / 6) * (np.sqrt(496 * x ** 2 + 232 * x + 1) - 22 * x + 1)
        # self.ax2.scatter(x, y)
        pareto_var_comb['y_analytic'] = y
        MSE = np.sum((pareto_var_comb['y_analytic'] - pareto_var_comb['y']) ** 2) / pareto_var_comb.shape[0]
        # print("MSE: ", MSE)
        self.ax3.scatter(n, MSE, c='r')

        ################################################################
        # calculate mean square error for averaged by grouping
        pareto_grp_avg = pareto_var_comb.groupby(
            pd.cut(df["x"], np.arange(self.x_bounds[0], self.x_bounds[1] + 0.01, 0.01))).mean()
        # get exact solution
        x = pareto_grp_avg['x']
        y = (1 / 6) * (np.sqrt(496 * x ** 2 + 232 * x + 1) - 22 * x + 1)
        # self.ax3.scatter(x, y)
        pareto_var_comb['y_analytic'] = y
        MSE_g = np.sum((pareto_grp_avg['y_analytic'] - pareto_grp_avg['y']) ** 2) / pareto_grp_avg.shape[0]
        self.ax3.scatter(n, MSE_g, facecolors="None", edgecolors='r', lw=2)

        if n + 1 == self.ng_max:
            self.ax1.scatter(pareto_shapes['f1'], pareto_shapes['f2'], marker='o', s=7, label=f'Pareto {n - 1}',
                             zorder=2)
            self.ax1.set_xlabel('f1')
            self.ax1.set_ylabel('f2')
            # self.ax2.scatter(pareto_var_comb['x'], pareto_var_comb['y'], marker='o', s=7, label=f'Pareto {n - 1}')
            self.ax2.set_xlabel('x')
            self.ax2.set_ylabel('y')
            # plot analytic solution
            xrange = np.linspace(*self.x_bounds, 50)
            yrange = np.linspace(*self.y_bounds, 50)
            X, Y = np.meshgrid(xrange, yrange)

            # F is one side of the equation, G is the other
            F = (8 * X + Y) * (6 * Y - 6)
            G = (2 * Y + X) * (2 * X - 2)

            self.ax2.contour(X, Y, (F - G), [0], linewidths=[3], colors=['k'])
            self.ax2.set_xlim(self.x_bounds)
            self.ax2.set_ylim(self.y_bounds)

            # plot analytic function
            f1 = self.f1(X, Y)
            f2 = self.f2(X, Y)

            self.ax1.scatter(f1, f2, s=5, color='y')

            # grouped data average
            self.ax2.scatter(pareto_grp_avg['x'], pareto_grp_avg['y'], marker='X', s=7, label=f'Pareto {n - 1}')
            self.ax3.set_yscale('log')
        # print(poc)

        df = df.loc[reorder_indx, :]
        self.df_global = df
        # self.recursive_save(df, "ga_res.xlsx", reorder_indx, poc)

        # crossover
        # print("Crossover")
        df_cross = self.crossover_sensitivity(df, n, self.fc)  # , elites["GR/Q
        # ic(df_cross)

        # mutation
        # print("Mutation")
        df_mutation = self.mutation(df, n, self.fm)
        # ic(df_mutation)

        # chaos
        # print("Chaos")
        df_chaos = self.chaos(self.fch, n)
        # ic(df_chaos)

        # take elites from previous generation over to next generation
        df_ng = pd.concat([df_cross, df_mutation, df_chaos], ignore_index=True)

        n += 1
        # print(n)
        # print("=" * 80)
        if n < self.ng_max:
            return self.run_sensitivity(n, df_ng)
        else:
            return

    @staticmethod
    def crossover(df, generation, f):  # , rq, grq
        elites = {}
        objectives = ['f1', 'f2']
        for o in objectives:
            elites[o] = df.sort_values(o)

        # naming convention G<generation number>_C<cavity number>_<type>
        # type refers to mutation M or crossover C
        df_co = pd.DataFrame(columns=["key", 'x', 'y'])

        for i in range(f):
            df_co.loc[i] = [f"G{generation}_C{i}_CO",
                            np.average(
                                [elites['f1'].loc[np.random.randint(30 if 30 < df.shape[0] else df.shape[0] - 1)]["x"],
                                 elites['f2'].loc[np.random.randint(30 if 30 < df.shape[0] else df.shape[0] - 1)][
                                     "x"]]),
                            # x
                            np.average(
                                [elites['f1'].loc[np.random.randint(30 if 30 < df.shape[0] else df.shape[0] - 1)]["y"],
                                 elites['f2'].loc[np.random.randint(30 if 30 < df.shape[0] else df.shape[0] - 1)][
                                     "y"]]),
                            ]
        return df_co

    @staticmethod
    def crossover_sensitivity(df, generation, f):  # , rq, grq
        elites = {}
        objectives = ['f1', 'f2']
        for o in objectives:
            elites[o] = df.sort_values(o)

        # naming convention G<generation number>_C<cavity number>_<type>
        # type refers to mutation M or crossover C
        df_co = pd.DataFrame(columns=["key", 'x', 'y'])

        for i in range(f):
            x_wavg, y_wavg = w_average(elites, [[0.895, 0.0973], [0.1, 0.9]], df.shape[0])
            df_co.loc[i] = [f"G{generation}_C{i}_CO",
                            x_wavg,  # x
                            y_wavg,
                            ]
        return df_co

    @staticmethod
    def mutation(df, n, f):
        # get range from table
        # get list based on mutation length
        if df.shape[0] < f:
            ml = np.arange(df.shape[0])
        else:
            ml = np.arange(f)

        df_ng_mut = pd.DataFrame(columns=['key', 'x', 'y'])

        df_ng_mut['x'] = df.loc[ml, "x"] * random.uniform(0.9, 1.1)
        df_ng_mut['y'] = df.loc[ml, "y"] * random.uniform(0.9, 1.1)

        key1, key2 = [], []
        for i in range(len(df_ng_mut)):
            key1.append(f"G{n}_C{i}_M")

        df_ng_mut['key'] = key1

        return df_ng_mut

    @staticmethod
    def chaos(f, n, xb, yb):
        data = {'key': [f"G{n}_C{i}_CH" for i in range(f)],
                'x': random.sample(list(np.linspace(*xb, f * 10)), f),
                'y': random.sample(list(np.linspace(*yb, f * 10)), f)}

        df = pd.DataFrame.from_dict(data)
        return df

    @staticmethod
    def pareto_front(df):

        datapoints = df.loc[:, ['f1', 'f2']] * (-1)
        pareto = oapackage.ParetoDoubleLong()
        # ic(datapoints)
        for ii in range(0, datapoints.shape[0]):
            w = oapackage.doubleVector(tuple(datapoints.iloc[ii].values))
            pareto.addvalue(w, ii)
        # pareto.show(verbose=1)  # Prints out the results from pareto

        lst = pareto.allindices()  # the indices of the Pareto optimal designs
        poc = len(lst)
        reorder_idx = list(lst) + [i for i in range(len(df)) if i not in lst]

        return reorder_idx, poc, lst

    def recursive_save(self, df, filename, pareto_index, poc):
        styler = self.color_pareto(df, poc)
        try:
            styler.to_excel(filename)
            # df.to_excel(filename)
        except PermissionError:
            filename = filename.split('.xlsx')[0]
            filename = fr'{filename}_1.xlsx'
            self.recursive_save(df, filename, pareto_index, poc)

    @staticmethod
    def color_pareto(df, no_pareto_optimal):
        def color(row):
            # if row.isnull().values.any():
            if row[0] in df['key'].tolist()[0:no_pareto_optimal]:
                return ['background-color: #6bbcd1'] * len(row)
            return [''] * len(row)

        # Save Styler Object for Later
        styler = df.style
        # Apply Styles (This can be chained or on separate lines)
        styler.apply(color, axis=1)
        # Export the styler to excel
        return styler

    @staticmethod
    def first_men(n, xb, yb):

        # create dataframe to apply constraints easily
        xrange_ = np.linspace(*xb, n)
        yrange_ = np.linspace(*yb, n)
        df__ = pd.DataFrame(np.array([xrange_, yrange_]).T, columns=['x', 'y'])
        # df__ = df__.loc[(df__['y'] + 9 * df__['x'] >= 6) & (-df__['y'] + 9 * df__['x'] >= 1)]
        xrange = df__['x']
        yrange = df__['y']

        data = {'key': [f"G0_C{i}_P" for i in range(len(xrange))],
                'x': xrange.tolist(),
                'y': yrange.tolist()}

        df = pd.DataFrame.from_dict(data)
        return df

    @staticmethod
    def f1(x, y):
        # convert explicitly to float
        x = x.astype(float)
        y = y.astype(float)

        # R. Duvigneau A. Zerbinati, J.A. D´esid´eri, Comparison between MGDA and PAES for multi-objective optimization, Tech. Rep. 7667 (INRIA Research Report No., 2011)
        # http://hal.inria.fr/inria-00605423.

        # https://en.wikipedia.org/wiki/Test_functions_for_optimization Constr-Ex problem

        f1 = 4 * x ** 2 + y ** 2 + x * y
        # f1 = x
        # f1 = x
        return f1

    @staticmethod
    def f2(x, y):
        # convert explicitly to float
        x = x.astype(float)
        y = y.astype(float)
        f2 = (x - 1) ** 2 + 3 * (y - 1) ** 2
        # f2 = (1 + y) / x
        # f2 = (1+y)*np.exp(-x/(1+y))
        return f2

    @staticmethod
    def show():
        # plt.legend()
        plt.show()


class GAMult:
    def __init__(self, n):
        self.pareto_shapes_sorted = None
        self.pareto_shapes = None
        self.u_bounds = None
        self.l_bounds = None
        self.var_names = None
        self.f2_interp = None
        self.df_global = pd.DataFrame()

        self.objs_dict = {}
        self.constraints_dict = {}

        self.ng_max = n
        self.objectives = []
        self.constraints = []
        self.weights = []
        self.processes = []
        self.df = None
        self.sbd = {}  # shape bounds dictionary
        self.poc = 5  # pareto optimal count

        # interpolation
        self.n_interp = 1000
        self.interp_error = []
        self.interp_error_avg = []

    def ea(self, vars_in, objs, constraints=None, func_constraints=None, n=0):
        """

        Parameters
        ----------
        vars_in: [['x1', [1, 2]], ['x3', [2, 4]], ..., ['xn', [1, 3.5]]]
        objs: [['equal', f1, '5'], ['min', f2], ['max', f3], ..., ['max', fn]]
        n: number of generations
        constraints: ['f1 < 20.0', 'f2 < 80.0', ..., 'fn > 20.0']

        Returns
        -------

        """
        if constraints is None:
            constraints = []

        if n == 0:
            # update lists
            # self.ng_max = 5
            self.objectives = objs
            self.vars_in = vars_in
            self.constraints = constraints
            self.func_constraints = func_constraints

            self.f2_interp = [np.zeros(self.n_interp) for _ in range(len(self.objectives))]

            # get input variable bounds
            self.var_names, self.l_bounds, self.u_bounds = [], [], []
            for v in vars_in:
                self.var_names.append(v[0])
                self.l_bounds.append(v[1][0])
                self.u_bounds.append(v[1][1])

                # include variable constraint to list of constraints
                self.constraints.append(f'{v[0]} >= {v[1][0]}')
                self.constraints.append(f'{v[0]} <= {v[1][1]}')

            self.df = self.generate_first_men(self.var_names, self.l_bounds, self.u_bounds, method='LHS')

        # evaluate objective
        for obj in objs:
            f = obj[1]
            self.df[f'{f.__name__}'] = f(*[self.df[x] for x in obj[2]])
            self.weights.append(1)

        # apply constraints
        if self.constraints:
            self.df = self.apply_bound_constraints(self.df, self.constraints)

        if self.func_constraints:
            self.df = self.apply_function_constraints(self.df, self.func_constraints)

        # optimize by page rank
        # remove the lowest ranking members
        df = self.df

        # compare with global dict and remove duplicates
        compared_cols = self.var_names
        if not self.df_global.empty:
            df = df.loc[~df.set_index(compared_cols).index.isin(
                self.df_global.set_index(compared_cols).index)]  # this line of code removes duplicates

        # naming convention G<generation number>_C<cavity number>_<type>
        # type refers to mutation M or crossover C

        # update with global dataframe
        if not self.df_global.empty:
            df = pd.concat([self.df_global, df], ignore_index=True)

            # # drop duplicates
            # df = df.drop_duplicates(subset=['A', 'B', 'a', 'b', 'Ri', 'L'], keep='first')

        # reset total rank
        df['total_rank'] = 0

        # rank shapes by objectives
        for i, obj in enumerate(self.objectives):
            if obj[0] == "min":
                df[f'rank_{obj[1].__name__}'] = df[obj[1].__name__].rank() * self.weights[i]
            elif obj[0] == "max":
                df[f'rank_{obj[1].__name__}'] = df[obj[1].__name__].rank(ascending=False) * self.weights[i]
            elif obj[0] == "equal":  # define properly later
                continue

            # if 'total_rank' in df.columns:
            df[f'total_rank'] = df[f'total_rank'] + df[f'rank_{obj[1].__name__}']

        # reorder
        tot = df.pop(f'total_rank')
        df[f'total_rank'] = tot / sum(self.weights)  # normalize by sum of weights

        # order shapes by rank
        df = df.sort_values(by=['total_rank'])
        df = df.reset_index(drop=True)

        # pareto condition
        reorder_indx, pareto_list = self.pareto_front(df)
        self.pareto_shapes = df.loc[pareto_list, [obj[1].__name__ for obj in self.objectives]]

        # estimate convergence
        obj_error = []
        obj0 = self.objectives[0][1]
        for i, obj in enumerate(self.objectives):
            if i != 0:
                # pareto_shapes = df.loc[:, [obj0.__name__, obj[1].__name__]]

                self.pareto_shapes_sorted = self.pareto_shapes.sort_values(obj0.__name__)
                f1 = np.linspace(min(self.pareto_shapes[obj0.__name__]), max(self.pareto_shapes[obj0.__name__]),
                                 self.n_interp)

                f2_interp = np.interp(f1, self.pareto_shapes_sorted[obj0.__name__],
                                      self.pareto_shapes_sorted[obj[1].__name__])
                rel_error = np.linalg.norm(f2_interp - self.f2_interp[i]) / max(np.abs(f2_interp))
                obj_error.append(rel_error)

                self.f2_interp[i] = f2_interp
        # ic(self.pareto_shapes_sorted)
        self.interp_error.append(max(obj_error))
        # ic(max(obj_error), np.average(self.interp_error))
        self.interp_error_avg.append(np.average(self.interp_error))
        # ic(self.interp_error, self.interp_error_avg)

        df = df.loc[reorder_indx, :]
        # reset index
        df = df.dropna().reset_index(drop=True)

        self.df_global = df

        # check if df_global is empty
        if self.df_global.shape[0] == 0:
            ic("Unfortunately, none survived the constraints and the program has to end. "
               "I can't even say that this was a good run.")
            return

        # save dataframe
        self.recursive_save(df, "../../utils/ga_res.xlsx", reorder_indx)

        # birth next generation
        # crossover
        # print("Crossover")
        df_cross = self.crossover(df, self.var_names, n, 100)

        # mutation
        # print("Mutation")
        df_mutation = self.mutation(df, self.var_names, n, 100)

        # chaos
        # print("Chaos")
        df_chaos = self.chaos(self.var_names, self.l_bounds, self.u_bounds, 'Uniform', 100, n)
        df_chaos_lhs = self.chaos(self.var_names, self.l_bounds, self.u_bounds, 'LHS', 100, n)
        df_chaos_rand = self.chaos(self.var_names, self.l_bounds, self.u_bounds, 'Random', 100, n)

        # take elites from previous generation over to next generation
        df_ng = pd.concat([df_cross, df_mutation, df_chaos, df_chaos_lhs, df_chaos_rand], ignore_index=True)

        # apply constraints
        if self.constraints:
            self.df = self.apply_bound_constraints(df_ng, self.constraints)

        # update dictionary
        self.df = df_ng

        n += 1
        ic(n)
        ic("=" * 80)
        if n < self.ng_max:
            return self.ea(vars_in, objs, n=n)
        else:
            self.n = n
            ic(df)
            return

    def apply_bound_constraints(self, df, constraints):
        # filter by constraints
        for const in constraints:
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

        return df

    def apply_function_constraints(self, df, constraints):
        # modify constraint
        constraints_mod = []
        for constraint in constraints:
            for obj in self.objectives:
                constraint = constraint.replace(obj[1].__name__, f'df["{obj[1].__name__}"]')
            constraints_mod.append(f'({constraint})')

        constraints_mod = ' & '.join(constraints_mod)
        # ic(constraints_mod)
        # apply constraint
        df = df.loc[eval(constraints_mod)]

        return df

    def weighted_mean_obj(self, tab_var, weights):
        rows_sims_no, cols = np.shape(tab_var)
        no_weights, dummy = np.shape(weights)  # z funckji quadr_stroud wekt columnowy

        if rows_sims_no == no_weights:
            expe = np.zeros((cols, 1))
            outvar = np.zeros((cols, 1))
            for i in range(cols):
                expe[i, 0] = np.dot(tab_var[:, i], weights)
                outvar[i, 0] = np.dot(tab_var[:, i] ** 2, weights)
            stdDev = np.sqrt(outvar - expe ** 2)
        else:
            expe = 0
            stdDev = 0
            ic('Cols_sims_no != No_weights')

        return list(expe.T[0]), list(stdDev.T[0])

    def generate_first_men(self, var_names, l_bounds, u_bounds, method='LHS', no_of_initial_points=300,
                           generation=0):
        columns = var_names
        dim = len(columns)

        if method == "LHS":
            const_var = []
            for i in range(dim - 1, -1, -1):
                if l_bounds[i] == u_bounds[i]:
                    const_var.append([columns[i], l_bounds[i]])
                    del columns[i]
                    l_bounds = np.delete(l_bounds, i)
                    u_bounds = np.delete(u_bounds, i)

            reduced_dim = len(columns)
            sampler = qmc.LatinHypercube(d=reduced_dim)
            _ = sampler.reset()
            sample = sampler.random(n=no_of_initial_points)
            # ic(qmc.discrepancy(sample))

            sample = qmc.scale(sample, l_bounds, u_bounds)

            df = pd.DataFrame()
            df['key'] = [f"G{generation}_C{i}_P" for i in range(no_of_initial_points)]
            df[columns] = sample

            for i in range(len(const_var) - 1, -1, -1):
                df[const_var[i][0]] = np.ones(no_of_initial_points) * const_var[i][1]

            return df

        # elif method == "Sobol Sequence":
        #     columns = ['A', 'B', 'a', 'b', 'Ri', 'L', 'Req']
        #     dim = len(columns)
        #     index = self.ui.sb_Sobol_Sequence_Index.value()
        #     l_bounds = np.array(list(self.sbd.values()))[:, 0]
        #     u_bounds = np.array(list(self.sbd.values()))[:, 1]
        #
        #     const_var = []
        #     for i in range(dim - 1, -1, -1):
        #         if l_bounds[i] == u_bounds[i]:
        #             const_var.append([columns[i], l_bounds[i]])
        #             del columns[i]
        #             l_bounds = np.delete(l_bounds, i)
        #             u_bounds = np.delete(u_bounds, i)
        #
        #     reduced_dim = len(columns)
        #     sampler = qmc.Sobol(d=reduced_dim, scramble=False)
        #     _ = sampler.reset()
        #     sample = sampler.random_base2(m=index)
        #     # ic(qmc.discrepancy(sample))
        #     # ic(l_bounds, u_bounds)
        #     sample = qmc.scale(sample, l_bounds, u_bounds)
        #
        #     df = pd.DataFrame()
        #     df['key'] = [f"G0_C{i}_P" for i in range(initial_points)]
        #     df[columns] = sample
        #
        #     for i in range(len(const_var) - 1, -1, -1):
        #         df[const_var[i][0]] = np.ones(initial_points) * const_var[i][1]
        #
        #     df['alpha_i'] = np.zeros(initial_points)
        #     df['alpha_o'] = np.zeros(initial_points)
        #
        #     return df
        elif method == "Random":
            X = list(zip(l_bounds, u_bounds))
            df = pd.DataFrame(np.array(
                [random.sample(list(np.linspace(*x, no_of_initial_points * 10)), no_of_initial_points) for x in X]).T,
                              columns=self.var_names)
            df['key'] = [f"G{generation}_C{i}_P" for i in range(no_of_initial_points)]
            return df
        elif method == "Uniform":
            # create dataframe to apply constraints easily
            X = list(zip(l_bounds, u_bounds))
            df = pd.DataFrame(np.array([np.linspace(*x, no_of_initial_points) for x in X]).T, columns=self.var_names)

            df['key'] = [f"G{generation}_C{i}_P" for i in range(no_of_initial_points)]
            return df

    def crossover(self, df, vars_names, generation, f):
        elites = {}
        for i, o in enumerate(self.objectives):
            if o[0] == "min":
                elites[f'{o[1].__name__}'] = df.sort_values(f'{o[1].__name__}')
            elif o[0] == "max":
                elites[f'{o[1].__name__}'] = df.sort_values(f'{o[1].__name__}', ascending=False)
            elif o[0] == "equal":
                continue

        obj_dict = {}
        for o in self.objectives:
            if o[0] != 'equal':
                obj_dict[o[1].__name__] = elites[o[1].__name__]

        obj = {}
        for key, o in obj_dict.items():
            obj[key] = o.reset_index(drop=True)

        # naming convention G<generation number>_C<cavity number>_<type>
        # type refers to mutation M or crossover C
        df_co = pd.DataFrame(columns=['key', *vars_names])

        for i in range(f):
            df_co.loc[i] = [f"G{generation}_C{i}_CO",
                            *[np.average([elites[f'{o[1].__name__}'].loc[
                                              np.random.randint(30 if 30 < df.shape[0] else df.shape[0] - 1)][f"{var}"]
                                          for o in self.objectives]) for var in vars_names]
                            ]
        return df_co

    def mutation(self, df, var_names, generation, f):
        # get range from table
        # get list based on mutation length
        if df.shape[0] < f:
            ml = np.arange(df.shape[0])
        else:
            ml = np.arange(f)

        df_ng_mut = pd.DataFrame(columns=['key', *var_names])

        for var in var_names:
            df_ng_mut[var] = df.loc[ml, var] * random.uniform(0.9, 1.1)

        key1, key2 = [], []
        for i in range(len(df_ng_mut)):
            key1.append(f"G{generation}_C{i}_M")

        df_ng_mut['key'] = key1

        return df_ng_mut

    def chaos(self, var_names, l_bounds, u_bounds, method='LHS', f=100, generation=0):
        df = self.generate_first_men(var_names, l_bounds, u_bounds, method=method, no_of_initial_points=f,
                                     generation=generation)
        return df

    @staticmethod
    def remove_duplicate_values(d):
        temp = []
        res = dict()
        for key, val in d.items():
            if val not in temp:
                temp.append(val)
                res[key] = val
        return res

    @staticmethod
    def proof_filename(filepath):
        # check if extension is included
        if filepath.split('.')[-1] != 'json':
            filepath = f'{filepath}.json'

        return filepath

    def recursive_save(self, df, filename, pareto_index):
        styler = self.color_pareto(df, self.poc)
        try:
            styler.to_excel(filename)
            # df.to_excel(filename)
        except PermissionError:
            filename = filename.split('.xlsx')[0]
            filename = fr'{filename}_1.xlsx'
            self.recursive_save(df, filename, pareto_index)

    def pareto_front(self, df):

        # datapoints = np.array([reverse_list(x), reverse_list(y), reverse_list(z)])
        # reverse list or not based on objective goal: minimize or maximize
        # datapoints = [self.negate_list(df.loc[:, o[1]], o[0]) for o in self.objectives]

        datapoints = df.loc[:, [o[1].__name__ for o in self.objectives]]
        # ic(datapoints)

        # ic(datapoints)
        for o in self.objectives:
            if o[0] == 'min':
                datapoints[o[1].__name__] = datapoints[o[1].__name__] * (-1)
            elif o[0] == "equal":
                pass

        # ic(datapoints)
        # convert datapoints to numpy array
        pareto = oapackage.ParetoDoubleLong()
        # ic(datapoints)

        for ii in range(0, datapoints.shape[0]):
            w = oapackage.doubleVector(tuple(datapoints.iloc[ii].values))
            pareto.addvalue(w, ii)
        pareto.show(verbose=1)  # Prints out the results from pareto

        lst = pareto.allindices()  # the indices of the Pareto optimal designs
        # ic(lst)
        self.poc = len(lst)
        # ic(self.poc)
        reorder_idx = list(lst) + [i for i in range(len(df)) if i not in lst]
        # ic(reorder_idx)
        # ic(lst)
        # optimal_datapoints = df.loc[lst, :]
        # ic(optimal_datapoints)
        # ic([optimal_datapoints.loc[i, :] for i in range(len(lst))])

        # return [optimal_datapoints[i, :] for i in range(datapoints.shape[0])]
        return reorder_idx, lst

    @staticmethod
    def negate_list(ll, arg):
        if arg == 'max':
            return ll  # to find the pareto maxima
        else:
            return [-x for x in ll]  # to find the pareto minima

    @staticmethod
    def overwriteFolder(invar, projectDir):
        path = fr"{projectDir}\SimulationData\SLANS\_process_{invar}"

        if os.path.exists(path):
            shutil.rmtree(path)
            dir_util._path_created = {}

        os.makedirs(path)

    @staticmethod
    def process_interval(interval_list):
        interval = []
        for i in range(len(interval_list) - 1):
            interval.append([interval_list[i], interval_list[i + 1]])

        return interval

    @staticmethod
    def color_pareto(df, no_pareto_optimal):
        def color(row):
            # if row.isnull().values.any():
            if row[0] in df['key'].tolist()[0:no_pareto_optimal]:
                return ['background-color: #6bbcd1'] * len(row)
            return [''] * len(row)

        # Save Styler Object for Later
        styler = df.style
        # Apply Styles (This can be chained or on separate lines)
        styler.apply(color, axis=1)
        # Export the styler to excel
        return styler

    def plot(self, title=''):

        # plot analytic solution
        xranges = {}
        for var in self.vars_in:
            xranges[var[0]] = np.linspace(*var[1], 50)

        X = np.meshgrid(*[xranges[var[0]] for var in self.vars_in])

        # plot analytic function
        functions = {}
        for obj in self.objectives:
            functions[obj[1].__name__] = obj[1](*[X[i] for i, _ in enumerate(obj[2])]).flatten()

        if len(self.objectives) == 2:
            fig, axs = plt.subplot_mosaic([[1, 3]])

            df_analytical = pd.DataFrame(functions)
            # apply function constraint
            if self.func_constraints:
                df_analytical = self.apply_function_constraints(df_analytical, self.func_constraints)

            axs[1].scatter(df_analytical['f1'], df_analytical['f2'], s=5, color='k', edgecolor='k', facecolor='none')

            axs[1].scatter(*[self.pareto_shapes_sorted[f"{obj[1].__name__}"] for obj in self.objectives], marker='o',
                           s=7, color='r',
                           label=f'Pareto Front (Gen. {self.n - 1})', zorder=2)
            axs[1].set_xlabel('$f_1$')
            axs[1].set_ylabel('$f_2$')
            axs[3].plot(self.interp_error, marker='P', label=r'$\epsilon$')
            # axs[3].plot(self.interp_error_avg, marker='X')
            axs[3].set_yscale('log')
            axs[1].legend()
            axs[3].legend()
            fig.suptitle(title)
            plt.show()
        elif len(self.objectives) == 3:
            fig = plt.figure()
            ax1 = fig.add_subplot(1, 1, 1, projection='3d')
            # ax2 = fig.add_subplot(2, 1, 2)

            df_analytical = pd.DataFrame(functions)
            # apply function constraint
            if self.func_constraints:
                df_analytical = self.apply_function_constraints(df_analytical, self.func_constraints)

            ax1.scatter(df_analytical['f1'], df_analytical['f2'], df_analytical['f3'], s=5, color='r', edgecolor='k', facecolor='none',)

            ax1.scatter(*[self.pareto_shapes_sorted[f"{obj[1].__name__}"] for obj in self.objectives], c='r', marker='o', s=20,
                           label=f'Pareto Front (Gen. {self.n - 1})', zorder=20)
            ax1.set_xlabel('$f_1$')
            ax1.set_ylabel('$f_2$')
            ax1.set_zlabel('$f_3$')
            ax1.set_xlim(0, 8)
            ax1.set_ylim(14.5, 17)
            ax1.set_zlim(-0.1, 0.2)
            ax1.view_init(elev=40., azim=-135)

            # ax2.plot(self.interp_error, marker='P')
            # ax2.plot(self.interp_error_avg, marker='X')
            # ax2.set_yscale('log')
            fig.suptitle(title)
            plt.legend(loc='upper right')
            plt.tight_layout()
            plt.show()


def w_average(df, weights, shape):
    x_wavg = weights[0][0] * df['f1'].loc[np.random.randint(30 if 30 < shape else shape - 1)]["x"] + \
             weights[1][0] * df['f2'].loc[np.random.randint(30 if 30 < shape else shape - 1)]["x"]

    y_wavg = weights[0][1] * df['f1'].loc[np.random.randint(30 if 30 < shape else shape - 1)]["y"] + \
             weights[1][1] * df['f2'].loc[np.random.randint(30 if 30 < shape else shape - 1)]["y"]

    return x_wavg, y_wavg


def f1(x, y, z):
    # convert explicitly to float
    x = x.astype(float)
    y = y.astype(float)

    # R. Duvigneau A. Zerbinati, J.A. D´esid´eri, Comparison between MGDA and PAES for multi-objective optimization, Tech. Rep. 7667 (INRIA Research Report No., 2011)
    # http://hal.inria.fr/inria-00605423.

    # https://en.wikipedia.org/wiki/Test_functions_for_optimization Constr-Ex problem

    # f = 4 * x ** 2 + y ** 2 + x * y
    # f = x
    # f = x
    f = 0.5 * (x ** 2 + y ** 2) + np.sin(x ** 2 + y ** 2)
    f = -10*(np.exp(-0.2*np.sqrt(x**2 + y**2)) + np.exp(-0.2*np.sqrt(y**2 + z**2)))
    return f


def f2(x, y, z):
    # convert explicitly to float
    x = x.astype(float)
    y = y.astype(float)
    # f = (x - 1) ** 2 + 3 * (y - 1) ** 2
    # f = (1 + y)/x
    # f = (1+y)*np.exp(-x/(1+y))
    f = (3 * x - 2 * y + 4)**2 / 8 + (x - y + 1)**2 / 27 + 15
    f = np.abs(x)**0.8 + 5*np.sin(x**3) + np.abs(y)**0.8 + 5*np.sin(y**3) + np.abs(z)**0.8 + 5*np.sin(z**3)
    return f


def f3(x, y):
    # convert explicitly to float
    x = x.astype(float)
    y = y.astype(float)
    f = 1 / (x ** 2 + y ** 2 + 1) - 1.1 * np.exp(-(x ** 2 + y ** 2))
    return f


if __name__ == '__main__':
    plot_settings()
    # ga = GA()
    # ga.run(1)
    # # ga.run_sensitivity(1)
    # ga.show()

    # ++++++++++++++++++++++++++++++++++++++++++
    # vars_in = [['x1', [-3, 3]], ['x2', [-3, 3]]]
    # objs = [['min', f1, ['x1', 'x2']], ['min', f2, ['x1', 'x2']], ['min', f3, ['x1', 'x2']]] # viennet
    # ++++++++++++++++++++++++++++++++++++++++++
    vars_in = [['x1', [-5, 5]], ['x2', [-5, 5]], ['x3', [-5, 5]]]
    objs = [['min', f1, ['x1', 'x2', 'x3']], ['min', f2, ['x1', 'x2', 'x3']]]
    n = 50
    constraints = None
    # func_constraints = ['f2/(0.858*np.exp(-0.541*f1)) >= 1', 'f2/(0.728*np.exp(-0.295*f1)) >= 1']
    # func_constraints = ['f2 + 9*f1 >= 6', '-f2 + 9*f1 >= 1']

    ga = GAMult(n)
    ga.ea(vars_in, objs, constraints)

    # ga.plot(title=r'$\mathbf{Viennet~function}$' + '\n' + r'Minimise = $\left\{ f_1(x, y) = 0.5\left(x^2+y^2\right)+\sin \left(x^2+y^2 \right) \right.,$ ' + '\n \t' +
    #               r'$f_2(x, y) = \frac{(3 x-2 y+4)^2}{8}+\frac{(x-y+1)^2}{27}+15,$' + '\n\t' +
    #               r'$\left. f_3(x, y) = \frac{1}{x^2+y^2+1}-1.1 \exp \left( -\left(x^2+y^2\right) \right) \right\}$')

    ga.plot(title=r'$\mathbf{Kursawe~function}$' + '\n' + r'Minimise = $\left\{ f_1(x_1, x_2, x_3) = \sum_{i=1}^2 \left[-0.2\sqrt{x_i^2+x_{i+1}^2}\right]\right.,$ ' + '\n \t' +
                  r'$\left. f_2(x_1, x_2, x_3) = \sum_{i=1}^3\left[\vert x_i \vert^{0.8} + 5\sin(x_i^3) \right] \right\}$')
