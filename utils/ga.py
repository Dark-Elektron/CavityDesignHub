import random
import numpy as np
import oapackage
import pandas as pd
from icecream import ic
from matplotlib import pyplot as plt


class GA:
    plt.rcParams["figure.figsize"] = (12, 6)

    def __init__(self):
        self.ng_max = 10
        self.df_global = None

        self.fig = plt.figure()

        gs = self.fig.add_gridspec(2, 2)
        self.ax1 = self.fig.add_subplot(gs[0, 0])
        self.ax2 = self.fig.add_subplot(gs[0, 1])
        # self.ax4 = self.fig.add_subplot(gs[1, 0])
        self.ax3 = self.fig.add_subplot(gs[1, :])

        self.x_bounds, self.y_bounds = [0.01, 1.5], [0.45, 1.5]
        # self.x_bounds, self.y_bounds = [0.1, 1], [0, 5]
        self.fm, self.fc, self.fch = 100, 100, 100

    def run(self, n, df_ng=None):
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

        ###############################################################
        # get exact solution
        x = pareto_var_comb['x']
        y = (1 / 6) * (np.sqrt(496 * x ** 2 + 232 * x + 1) - 22 * x + 1)
        # self.ax3.scatter(x, y)
        pareto_var_comb['y_analytic'] = y
        MSE = np.sum((pareto_var_comb['y_analytic'] - pareto_var_comb['y']) ** 2) / pareto_var_comb.shape[0]
        # print("MSE: ", MSE)
        self.ax3.scatter(n, MSE, c='k')

        ################################################################
        # calculate mean square error for averaged by grouping
        pareto_grp_avg = pareto_var_comb.groupby(pd.cut(df["x"], np.arange(self.x_bounds[0], self.x_bounds[1] + 0.01, 0.01))).mean()
        # get exact solution
        x = pareto_grp_avg['x']
        y = (1 / 6) * (np.sqrt(496 * x ** 2 + 232 * x + 1) - 22 * x + 1)
        self.ax2.scatter(x, y)
        pareto_var_comb['y_analytic'] = y
        MSE_g = np.sum((pareto_grp_avg['y_analytic'] - pareto_grp_avg['y']) ** 2) / pareto_grp_avg.shape[0]
        self.ax3.scatter(n, MSE_g, facecolors="None", edgecolors='k', lw=2)

        if n + 1 == self.ng_max:
            self.ax1.scatter(pareto_shapes['f1'], pareto_shapes['f2'], marker='o', s=7, label=f'Pareto {n - 1}', zorder=2)
            self.ax1.set_xlabel('f1')
            self.ax1.set_ylabel('f2')
            # self.ax2.scatter(pareto_var_comb['x'], pareto_var_comb['y'], marker='o', s=7, label=f'Pareto {n - 1}')
            self.ax2.set_xlabel('x')
            self.ax2.set_ylabel('y')

            # plot analytic solution
            xrange = np.linspace(*self.x_bounds, 50)
            yrange = np.linspace(*self.y_bounds, 50)
            X, Y = np.meshgrid(xrange, yrange)
            idx = (Y + 9*X >= 6) & (-Y + 9*X >= 1)
            X[~idx] = 0
            Y[~idx] = 0

            # F is one side of the equation, G is the other
            F = (8 * X + Y) * (6 * Y - 6)
            G = (2 * Y + X) * (2 * X - 2)

            # self.ax2.contour(X, Y, (F - G), [0], linewidths=[3], colors=['k'])
            self.ax2.set_xlim(self.x_bounds)
            self.ax2.set_ylim(self.y_bounds)

            # plot analytic function
            f1 = self.f1(X, Y)
            f2 = self.f2(X, Y)

            self.ax1.scatter(f1, f2, s=5, color='r')

            # grouped data average
            self.ax2.scatter(pareto_grp_avg['x'], pareto_grp_avg['y'], marker='P', s=7, label=f'Pareto {n - 1}')
            self.ax3.set_yscale('log')

        # print(poc)

        df = df.loc[reorder_indx, :]
        self.df_global = df
        self.recursive_save(df, "ga_res.xlsx", reorder_indx, poc)

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
        df_chaos = self.chaos(self.fch, n)
        # ic(df_chaos)

        # take elites from previous generation over to next generation
        df_ng = pd.concat([df_cross, df_mutation, df_chaos], ignore_index=True)

        # apply constraints
        df_ng = df_ng.loc[(df_ng['x'] <= self.x_bounds[1]) & (df_ng['x'] >= self.x_bounds[0])]
        df_ng = df_ng.loc[(df_ng['y'] <= self.y_bounds[1]) & (df_ng['y'] >= self.y_bounds[0])]
        df_ng = df_ng.loc[(df_ng['y'] + 9*df_ng['x'] >= 6) & (-df_ng['y'] + 9*df_ng['x'] >= 1)]

        n += 1
        # print(n)
        # print("=" * 80)
        if n < self.ng_max:
            return self.run(n, df_ng)
        else:
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
        pareto_grp_avg = pareto_var_comb.groupby(pd.cut(df["x"], np.arange(self.x_bounds[0], self.x_bounds[1] + 0.01, 0.01))).mean()
        # get exact solution
        x = pareto_grp_avg['x']
        y = (1 / 6) * (np.sqrt(496 * x ** 2 + 232 * x + 1) - 22 * x + 1)
        # self.ax3.scatter(x, y)
        pareto_var_comb['y_analytic'] = y
        MSE_g = np.sum((pareto_grp_avg['y_analytic'] - pareto_grp_avg['y']) ** 2) / pareto_grp_avg.shape[0]
        self.ax3.scatter(n, MSE_g, facecolors="None", edgecolors='r', lw=2)

        if n + 1 == self.ng_max:
            self.ax1.scatter(pareto_shapes['f1'], pareto_shapes['f2'], marker='o', s=7, label=f'Pareto {n - 1}', zorder=2)
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
                                 elites['f2'].loc[np.random.randint(30 if 30 < df.shape[0] else df.shape[0] - 1)]["x"]]),
                            # x
                            np.average(
                                [elites['f1'].loc[np.random.randint(30 if 30 < df.shape[0] else df.shape[0] - 1)]["y"],
                                 elites['f2'].loc[np.random.randint(30 if 30 < df.shape[0] else df.shape[0] - 1)]["y"]]),
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

        df_ng_mut.loc[:, 'x'] = df.loc[ml, "x"] * random.uniform(0.9, 1.1)
        df_ng_mut.loc[:, 'y'] = df.loc[ml, "y"] * random.uniform(0.9, 1.1)

        key1, key2 = [], []
        for i in range(len(df_ng_mut)):
            key1.append(f"G{n}_C{i}_M")

        df_ng_mut.loc[:, 'key'] = key1

        return df_ng_mut

    @staticmethod
    def chaos(f, n):
        data = {'key': [f"G{n}_C{i}_CH" for i in range(f)],
                'x': random.sample(list(np.linspace(0.01, 1.5, f * 10)), f),
                'y': random.sample(list(np.linspace(0.45, 1.5, f * 10)), f)}

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
        pareto.show(verbose=1)  # Prints out the results from pareto

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
        df__ = df__.loc[(df__['y'] + 9 * df__['x'] >= 6) & (-df__['y'] + 9 * df__['x'] >= 1)]
        xrange = df__['x']
        yrange = df__['y']

        data = {'key': [f"G0_C{i}_P" for i in range(len(xrange))],
                'x': xrange.tolist(),
                'y': yrange.tolist()}

        return pd.DataFrame.from_dict(data)

    @staticmethod
    def f1(x, y):
        f1 = 4 * x ** 2 + y ** 2 + x * y
        # f1 = x
        return f1

    @staticmethod
    def f2(x, y):
        f2 = (x - 1) ** 2 + 3 * (y - 1) ** 2
        # f2 = (1 + y)/x
        return f2

    @staticmethod
    def show():
        plt.legend()
        plt.show()


def w_average(df, weights, shape):
    x_wavg = weights[0][0] * df['f1'].loc[np.random.randint(30 if 30 < shape else shape - 1)]["x"] + \
             weights[1][0] * df['f2'].loc[np.random.randint(30 if 30 < shape else shape - 1)]["x"]

    y_wavg = weights[0][1] * df['f1'].loc[np.random.randint(30 if 30 < shape else shape - 1)]["y"] + \
             weights[1][1] * df['f2'].loc[np.random.randint(30 if 30 < shape else shape - 1)]["y"]

    return x_wavg, y_wavg


if __name__ == '__main__':
    ga = GA()
    ga.run(1)
    # ga.run_sensitivity(1)
    ga.show()
