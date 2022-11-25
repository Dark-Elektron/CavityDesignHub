import itertools
from itertools import combinations
import matplotlib.ticker as mticker
from matplotlib import pyplot as plt
from scipy.stats import qmc
import numpy as np
import pandas as pd
from icecream import ic
from sklearn import linear_model
from sympy import symbols


class PCE:
    def __init__(self, df, rand_vars, obj_vars, poly_type):
        self.df = df
        self.rand_vars = rand_vars
        self.obj_vars = obj_vars
        poly_type_dict = {"legendre": "Le",
                          "hermite": "He",
                          "laguerre": "La",
                          "jacobi": "Ja"
                          }
        self.poly_type = poly_type_dict[poly_type.lower()]
        self.pce = None

    def get_pce(self, p_order=1, truncation=1):

        random_var_order = [p_order for v in self.rand_vars]
        # random_var_order = [1, 1, 1, 1]
        rvo = [[i for i in range(x+1)] for x in random_var_order]

        alpha = list(itertools.product(*rvo))
        x = {}
        for r in self.rand_vars:
            x[r] = self.df[r]

        poly_list, poly_chaos_ex, poly_list_sym, symbols_dict = {}, {}, {}, {}
        ic(len(alpha))

        for a in alpha:
            # truncate
            if sum(a) <= truncation:
                # if np.linalg.norm(a) <= truncation:
                ll = []
                for j, i in enumerate(a):
                    if self.poly_type.lower() == 'le':
                        poly_list[f'{a}'] = [self.Le(i, x[list(x.keys())[j]]) for j, i in enumerate(a)]
                        pc, sy = self.Le(i, x[list(x.keys())[j]], self.rand_vars[j])
                    elif self.poly_type.lower() == 'he':
                        poly_list[f'{a}'] = [self.He(i, x[list(x.keys())[j]]) for j, i in enumerate(a)]
                        pc, sy = self.He(i, x[list(x.keys())[j]], self.rand_vars[j])
                    elif self.poly_type.lower() == 'la':
                        poly_list[f'{a}'] = [self.Le(i, x[list(x.keys())[j]]) for j, i in enumerate(a)]
                        pc, sy = self.Le(i, x[list(x.keys())[j]], self.rand_vars[j])
                    elif self.poly_type.lower() == 'ja':
                        poly_list[f'{a}'] = [self.Le(i, x[list(x.keys())[j]]) for j, i in enumerate(a)]
                        pc, sy = self.Le(i, x[list(x.keys())[j]], self.rand_vars[j])
                    else:
                        print("The polynomial type input is not recognised. It is set to Legendre by default.")
                        poly_list[f'{a}'] = [self.Le(i, x[list(x.keys())[j]]) for j, i in enumerate(a)]
                        pc, sy = self.Le(i, x[list(x.keys())[j]], self.rand_vars[j])

                    ll.append(pc)

                    if len(symbols_dict) <= len(a):
                        symbols_dict[sy] = 0

                poly_list_sym[f'{a}'] = np.prod(ll)

                poly_chaos_ex[f'{a}'] = [f"pce_model.{self.poly_type}({i}, DF['{list(x.keys())[j]}'])"
                                         for j, i in enumerate(a)]

        # regeression
        _, _ = self.regression(poly_list_sym, symbols_dict)

        c = {}
        poly_sym_d = {}

        for ob in self.obj_vars:
            c[ob] = {}
            poly_sym = 0
            ci = 0
            for k, v in poly_list.items():
                p = np.prod(np.vstack(v), axis=0)

                c[ob][k] = np.dot(self.df[ob], p)/np.dot(p, p)
                poly_sym += c[ob][k]*poly_list_sym[k]

                ci += 1

            poly_sym_d[ob] = poly_sym

        # build polynomials
        obj_pce = {}
        for ob in self.obj_vars:
            obj_pce[ob] = ' + '.join(map(str, [f"{coeff}*{'*'.join(map(str, poly_chaos_ex[key]))}"
                                               for key, coeff in c[ob].items()]))

        self.pce = obj_pce

        return obj_pce, c

    def projection(self):
        return

    def regression(self, poly_list_sym, symbols_dict):
        reg = linear_model.LinearRegression(fit_intercept=True)
        # build matrix
        A = []
        # print(poly_list_sym)
        for k, v in self.df.iterrows():
            symbols_dict = self.update_symbols_dict(symbols_dict, v[self.rand_vars])
            aa = [p.subs(symbols_dict) if not isinstance(p, np.int32) else 1 for p in poly_list_sym.values()]
            A.append(aa)

        poly, coef_dict = {}, {}
        for ob in self.obj_vars:
            b = self.df[ob].to_numpy(dtype='float')
            reg.fit(np.array(A, dtype='float'), b.reshape(-1, 1))

            coef_dict[ob] = reg.coef_.copy()[0]
            coef_dict[ob][0] = reg.intercept_[0]
            ic(reg.score(A, b))
            # ic(x)

            poly[ob] = np.sum(np.array(list(poly_list_sym.values()))*coef_dict[ob])

            symbols_dict = self.update_symbols_dict(symbols_dict, self.df.loc[1, self.rand_vars])
            # ic(poly[ob].subs(symbols_dict))

        self.pce = poly
        ic(poly)

        return poly, coef_dict

    @staticmethod
    def Le(n, x, symbol=None):
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
            shift = (max(x) + min(x)) / 2
            scale = (max(x) - min(x)) / 2

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

    @staticmethod
    def He(n, x, symbol=None):
        shift = (max(x) + min(x)) / 2
        scale = (max(x) - min(x)) / 2
        if n == 0:
            return np.ones(len(x))

        if n == 1:
            x = (np.array(x) - shift) / scale
            return 2 * x

        if n == 2:
            x = (np.array(x) - shift) / scale
            return 4 * x ** 2 - 2

        if n == 3:
            x = (np.array(x) - shift) / scale
            return 8 * x ** 3 - 12 * x

        if n == 4:
            x = (np.array(x) - shift) / scale
            return 16 * x ** 4 - 48 * x ** 2 + 12

    def self_validation(self):
        # find least squares error
        pce_model = self
        DF = pd.DataFrame()  # it is also important that the name is DF # change later
        for key, vec in self.df.items():
            DF[key] = vec.to_numpy()

        # evaluate polynomial expansion on uq_model input dataframe
        if self.pce is not None:
            print(self.pce)
            for key, ob in self.pce.items():
                DF[key] = eval(ob)

        # find difference
        print(self.obj_vars)
        for obj in self.obj_vars:
            diff = self.df[obj] - DF[obj]
            # print(diff)
            print(f"error norm {obj}: ", np.linalg.norm(diff))

    def cross_validation(self):
        pass

    @staticmethod
    def update_symbols_dict(sym_dict, values):
        assert len(sym_dict) == len(values)
        new_sym_list = {}
        for i, a in enumerate(sym_dict.keys()):
            new_sym_list[a] = values[i]

        return new_sym_list


class UQModel:

    def __init__(self):
        self.input_variable_names = None
        self.objective_variables = None
        self.input_variable_bounds = None
        self.input_variable_distribution = None
        self.var_sample_size = 1000
        self.model = None
        self.model_input = None
        self.model_output = None

    def run_analysis(self):
        """Model is a function or routine that takes in some input variables
        as a pandas dataframe"""

        # df = self.generate_input_space()
        # self.model_output = self.model(self.model_input, df)

        df_ = self.generate_input_qmc(self.var_sample_size, self.input_variable_names, self.input_variable_bounds)
        model_output_satelli = self.model(self.model_input, df_)
        # Sj = self.sobol_satelli(model_output_satelli, self.input_variable_names, self.objective_variables, 10000)
        # self.sobol_janon(model_output_satelli, self.input_variable_names, self.objective_variables, 10000)

        self.model_output = model_output_satelli

        return self.model_output

    def set_model(self, model):
        self.model = model

    def method(self):
        pass

    def set_input_variables(self, dd):
        """Structure of dd
        dd = {
            "random variables": 'x1', 'x2', ..., 'xn',
            "objective variables": 'v1', 'v2', ..., 'vm',
            'bounds': [[]1, []2,..., []n],
            'distribution': [d1, d2, ..., dn]
            }
        }
        If the distribution is uniform, bounds holds a list of the lower and upper bounds of the corresponding
        variable, respectively.
        if the distribution is set to Lognormal, bounds holds a list of the mean and standard deviation of the
        corresponding variable, respectively
        """

        self.input_variable_names = dd['random variables']
        self.objective_variables = dd['objective variables']
        self.input_variable_bounds = dd['bounds']
        self.input_variable_distribution = dd['distribution']
        # check input variables that all input variable dimensions are equal
        assert len(self.input_variable_names) == len(self.input_variable_distribution)
        assert len(self.input_variable_names) == len(self.input_variable_bounds)

        self.model_input = dd

    def generate_input_space(self):
        var = self.input_variable_names
        bounds = np.array(self.input_variable_bounds)

        dd = {}

        for k, v in enumerate(var):
            if self.input_variable_distribution[k] == 'Uniform':
                dd[v] = np.linspace(min(bounds[k]), max(bounds[k]), int(self.var_sample_size))
            elif self.input_variable_distribution[k] == 'Lognormal':
                eta = np.sqrt(np.log((bounds[k][1]/bounds[k][0])**2 + 1))
                mean = np.log(bounds[k][0]) - eta**2/2
                # print(mean, eta)
                dd[v] = np.linspace(min(bounds[k]), max(bounds[k]), int(self.var_sample_size))
                # dd[v] = np.random.lognormal(bounds[k][0], bounds[k][1], self.var_sample_size)

        mg = np.meshgrid(*dd.values())

        df = pd.DataFrame()
        for i, rand_var in enumerate(var):
            df[rand_var] = mg[i].flatten()

        print(df.shape)

        return df

    @staticmethod
    def generate_input_mc(N, rand_vars, bounds):
        D = len(rand_vars)

        A = np.random.uniform(0, 1, (N, D))
        B = np.random.uniform(0, 1, (N, D))

        # scale
        for i in range(D):
            lb, ub = bounds[i]
            A[:, i] = lb + A[:, i]*(ub - lb)
            B[:, i] = lb + B[:, i]*(ub - lb)

        df = pd.DataFrame(A, columns=rand_vars)
        df = pd.concat([df, pd.DataFrame(B, columns=rand_vars)], ignore_index=True)

        for i in range(D):
            A_copy = np.copy(A)
            A_copy[:, i] = B[:, i]
            df = pd.concat([df, pd.DataFrame(A_copy, columns=rand_vars)], ignore_index=True)
        print(df.shape)
        return df

    @staticmethod
    def generate_input_qmc(N, rand_vars, bounds):
        D = len(rand_vars)
        l_bounds, u_bounds = np.array(bounds).T
        sampler = qmc.LatinHypercube(d=D)
        _ = sampler.reset()
        # generate matrix A
        A = sampler.random(n=N)
        A = qmc.scale(A, l_bounds, u_bounds)
        # generate matrix B
        B = sampler.random(n=N)
        B = qmc.scale(B, l_bounds, u_bounds)

        df = pd.DataFrame(A, columns=rand_vars)
        df = pd.concat([df, pd.DataFrame(B, columns=rand_vars)], ignore_index=True)

        for i in range(D):
            A_copy = np.copy(A)
            A_copy[:, i] = B[:, i]

            df = pd.concat([df, pd.DataFrame(A_copy, columns=rand_vars)], ignore_index=True)

        return df

    def set_sample_size(self, n):
        self.var_sample_size = n

    def example_simplySupportedBeam(self):
        self.input_variable_names = ['b', 'h', 'L', 'E', 'p']
        self.input_variable_bounds = np.array([[0.15, 0.0075], [0.3, 0.015], [5, 0.05], [3e10, 4.5e9], [1e4, 2e3]])
        n = self.var_sample_size

        var = self.input_variable_names
        var_int = self.input_variable_bounds
        # df = generate_sobol_sequence(5, 8, var, var_int)

        dd = {}
        for k, v in enumerate(var):
            dd[v] = np.linspace(min(var_int[k]), max(var_int[k]), int(n))
        mg = np.meshgrid(*dd.values())
        df = pd.DataFrame()
        for i, rand_var in enumerate(var):
            df[rand_var] = mg[i].flatten()

        b, h, L, E, p = df.to_numpy().T

        # The beam is considered primatic, therefore:
        Im = b * h ** 3 / 12  # the moment of inertia

        # now for the actual execution we use a vectorized formula:
        Y = np.zeros((len(L), 9))
        for jj in range(1, 10):
            # calculate the xi values:
            xi = jj / 10 * L

            Y[:, jj - 1] = -p * xi * (L ** 3 - 2 * xi ** 2 * L + xi ** 3) / (24 * E * Im)

        Y = pd.concat([df, pd.DataFrame(Y, columns=['Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8', 'Y9'])], axis=1)
        return Y

    @staticmethod
    def generate_sobol_sequence(dim, index, columns, bounds):
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

    @staticmethod
    def sub_lists(ll):
        # initializing empty list
        comb = []

        # Iterating till length of list
        for i in range(len(ll) + 1):
            # Generating sub list
            comb += [list(j) for j in combinations(ll, i)]
        # Returning list
        return comb[1:-1]

    def sobol_df(self):
        x = self.model_input['random variables']
        obj = self.model_input['objective variables']
        df_ = self.model_output

        power_set = self.sub_lists(x)

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
                    var_i = 1 - np.var(df_.groupby(subset_inv)[ob].mean()) / np.var(df_[ob])
                    STi[ob][f'V[E[{ob}|{",".join(map(str, subset))}]]'] = var_i
                else:
                    power_subset = self.sub_lists(subset)
                    for rand_var in power_subset:
                        var = var - Sj[ob][f'V[E[{ob}|{",".join(map(str, rand_var))}]]']

                    Sj[ob][f'V[E[{ob}|{",".join(map(str, subset))}]]'] = var

        return Sj, STi

    def sobol_satelli(self):
        """
            A starts from 0 to N-1
            B starts from N to 2N - 1
            AB0 starts from 2N to 3N - 1
        """
        rand_vars = self.model_input['random variables']
        obj_vars = self.model_input['objective variables']
        df = self.model_output
        N = self.var_sample_size

        Sj = {}
        STj = {}
        for ob in obj_vars:
            Sj[ob] = {}
            STj[ob] = {}
            for j, rv in enumerate(rand_vars):
                E2Y = (1/N)*(np.sum([df[ob][i] * df[ob][N + i] for i in range(N)]))
                VEx = (1/N)*np.sum([df[ob][N + i] * df[ob][(j+2)*N + i] for i in range(N)]) - E2Y

                EY = (1/N)*np.sum([(df[ob][i]) for i in range(N)])
                VY = (1/N)*np.sum([(df[ob][i])**2 for i in range(N)]) - EY**2

                EVYx = (1/N)*np.sum([(df[ob][i])**2 - df[ob][i]*df[ob][(j+2)*N + i] for i in range(N)])
                # EVYx = (1/(2*N))*np.sum([(df[ob][i] - df[ob][(j+2)*N + i])**2 for i in range(N)])

                Sj[ob][f"V[E[{ob}|{rv}]]"] = VEx/VY
                STj[ob][f"V[E[{ob}|{rv}]]"] = EVYx/VY

        ic(Sj)
        return Sj, STj

    def sobol_janon(self):
        """
            A starts from 0 to N-1
            B starts from N to 2N - 1
            AB0 starts from 2N to 3N - 1
        """
        rand_vars = self.model_input['random variables']
        obj_vars = self.model_input['objective variables']
        df = self.model_output
        N = self.var_sample_size

        Sj, STj = {}, {}

        for ob in obj_vars:
            Sj[ob], STj[ob] = {}, {}
            for j, rv in enumerate(rand_vars):
                E2Y = ((1/N)*(np.sum([(df[ob][N+i] + df[ob][(j+2)*N + i])/2 for i in range(N)])))**2
                VEx = (1/N)*np.sum([df[ob][N + i] * df[ob][(j+2)*N + i] for i in range(N)]) - E2Y

                EY = ((1/N)*(np.sum([(df[ob][N+i] + df[ob][(j+2)*N + i])/2 for i in range(N)])))
                VY = (1/N)*np.sum([(df[ob][N+i]**2 + df[ob][(j+2)*N + i]**2)/2 for i in range(N)]) - EY**2

                EVYx = (1/N)*np.sum([(df[ob][i])**2 - df[ob][i]*df[ob][(j+2)*N + i] for i in range(N)])
                # EVYx = (1/(2*N))*np.sum([(df[ob][i] - df[ob][(j+2)*N + i])**2 for i in range(N)])

                Sj[ob][f"V[E[{ob}|{rv}]]"] = VEx/VY
                STj[ob][f"V[E[{ob}|{rv}]]"] = EVYx/VY

        ic(Sj)
        return Sj, STj

    def plot_sobol(self, S, ylabel='Sobol', table=False, plot_type="Stacked"):
        # plot
        fig, ax = plt.subplots()
        x = self.input_variable_names
        obj = self.objective_variables
        width = 0.5

        # for ob in obj:
        # rand_var_dict = {0: 'A', 1: 'B', 2: 'a', 3: 'b', 4: 'Ri'}
        # rand_var_dict = {0: 'lh1', 1: 'lh3', 2: 'lh4', 3: 'dh3', 4: 'alpha_h', 5: 'ch2', 6: 'r_cyl', 7: 'offset_y'}
        rand_var_dict = {}
        for i, v in enumerate(x):
            rand_var_dict[i] = v

        if plot_type == "Stacked":
            bottom, bottom1 = np.zeros(len(obj)), np.zeros(len(obj))
            for key in x:
                # for key in x.keys():
                if key == 0:
                    ax.bar(obj, [S[ob][f'V[E[{ob}|{key}]]'] for ob in obj], width, label=r"$\mathbf{" + key + '}$')
                else:
                    ax.bar(obj, [S[ob][f'V[E[{ob}|{key}]]'] for ob in obj], width, label=r"$\mathbf{" + key + '}$',
                           bottom=bottom)

                bottom += np.array([S[ob][f'V[E[{ob}|{key}]]'] for ob in obj])
                # bottom1 += np.array([STi[ob][f'V[E[{ob}|{key}]]'] for ob in obj])

            ticks_loc = ax.get_xticks()
            ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            ax.set_xticklabels(obj)
            ax.set_ylabel(ylabel)

            if table:
                self.plot_table(S, x, ax)

            lines, labels = ax.get_legend_handles_labels()

            if int(len(rand_var_dict.keys())) > 5:
                ncol = 5
            else:
                ncol = int(len(rand_var_dict.keys()))
            fig.legend(lines, labels, loc="upper center",
                       ncol=ncol,
                       fancybox=True, shadow=False,
                       bbox_to_anchor=(0.5, 1.0))

            # ax[0].legend(loc="upper center", ncol=int(len(rand_var_dict.keys())), bbox_to_anchor=(0.0, 1.1),
            #        fancybox=True, shadow=True)
            # ax[1].legend(loc="upper center", ncol=int(len(rand_var_dict.keys())), bbox_to_anchor=(0.0, 1.1),
            #           fancybox=True, shadow=True)
        else:
            data = []
            X = np.arange(len(x))
            width = 1/(len(x) + 10)
            for i, ob in enumerate(obj):
                data = list(S[ob].values())
                ax.bar(X + i*width, data, width=width, label=ob)

            plt.xticks([r + width*len(x) for r in range(len(x))], x)

            lines, labels = ax.get_legend_handles_labels()
            if len(x) > 5:
                ncol = 5
            else:
                ncol = len(x)
            fig.legend(lines, labels, loc="upper center",
                       ncol=ncol,
                       fancybox=True, shadow=False,
                       bbox_to_anchor=(0.5, 1.0))

        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_table(S, x, ax):
        # Add a table at the bottom of the axes
        cellText = []
        rowLabels = []
        colLabels = []
        temp_col = []
        n = 0
        for k, v in S.items():
            rowLabels.append(r"$\mathbf{" + k + '}$')

            temp_col, temp_cellText = [], []
            for xx, xv in zip(list(v.keys()), list(v.values())):
                if xx in [f'V[E[{k}|{p}]]' for p in x]:
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
        the_table = ax.table(cellText=cellText,
                             rowLabels=rowLabels,
                             colLabels=colLabels,
                             rowColours=rcolors,
                             colColours=ccolors,
                             # cellColours=colours,
                             loc='bottom', bbox=[0.0, -0.35, 1, .28])
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(8)

        plt.subplots_adjust(bottom=0.3)
        plt.subplots_adjust(top=0.9)


def func(model_input, df):
    x = model_input['random variables']
    obj = model_input['objective variables']

    # The beam is considered primatic, therefore:
    Y = np.zeros((df.shape[0], len(obj)))

    for jj, ob in enumerate(obj):
        Y[:, jj-1] = (1 - df[x[0]]) ** 2 + 100 * (df[x[1]] - df[x[0]] ** 2) ** 2

    Y = pd.concat([df, pd.DataFrame(Y, columns=obj)], axis=1)

    return Y


def ishigami(model_input, df):
    x = model_input['random variables']
    obj = model_input['objective variables']
    # The beam is considered primatic, therefore:
    Y = np.zeros((df.shape[0], len(obj)))

    for jj, ob in enumerate(obj):
        a, b = 7, 0.1
        Y[:, jj - 1] = np.sin(df[x[0]]) + a*(np.sin(df[x[1]]))**2 + b*df[x[2]]**4*np.sin(df[x[0]])

    Y = pd.concat([df, pd.DataFrame(Y, columns=obj)], axis=1)

    return Y


def example_simplySupportedBeam(model_input, df):

    x = model_input['random variables']
    obj = model_input['objective variables']

    b, h, L, E, p = df[x[0]], df[x[1]], df[x[2]], df[x[3]], df[x[4]]

    # The beam is considered primatic, therefore:
    Im = b * h ** 3 / 12  # the moment of inertia

    # now for the actual execution we use a vectorized formula:
    Y = np.zeros((df.shape[0], len(obj)))
    for jj in range(1, 10):
        # calculate the xi values:
        xi = jj / 10 * L

        Y[:, jj - 1] = -p * xi * (L ** 3 - 2 * xi ** 2 * L + xi ** 3) / (24 * E * Im)

    Y = pd.concat([df, pd.DataFrame(Y, columns=obj)], axis=1)
    # print(Y)
    return Y


def toy_problem(model_input, df):
    x = model_input['random variables']
    obj = model_input['objective variables']

    x, y = df[x[0]], df[x[1]]

    f1 = 4 * x ** 2 + y ** 2 + x * y
    f2 = (x-1)**2 + 3*(y-1)**2
    F = np.array([f1, f2]).T

    F = pd.concat([df, pd.DataFrame(F, columns=obj)], axis=1)

    return F


def dqw_pce_from_df(model_input, df):
    # create function from data input
    filename = fr'D:\CST Studio\Hook Coupler Study\3. Optimize Hook Coupler Geometry\DQW_Smax_Fmax_Data_11var_2p.xlsx'
    df_data = pd.read_excel(filename, 'Sheet1')
    x = model_input['random variables']
    obj = model_input['objective variables']

    pce_model = PCE(df_data, x, obj, 'Legendre')  # it is important that the name is pce_model
    pce, pce_coeff = pce_model.get_pce(2, 2)
    pce_model.self_validation()

    DF = pd.DataFrame()  # it is also important that the name is DF # change later
    for key, vec in df.items():
        DF[key] = vec.to_numpy()

    # evaluate polynomial expansion on uq_model input dataframe
    for key, ob in pce.items():
        DF[key] = eval(ob)

    return DF


if __name__ == '__main__':

    dd = {
        "random variables": ['x1', 'x2'],
        "objective variables": ['v1'],
        'bounds': [[-2, 2], [-2, 2]],
        'distribution': ['Uniform', 'Uniform']
    }
    dd_ishigami = {
        "random variables": ['x1', 'x2', 'x3'],
        "objective variables": ['v1'],
        'bounds': [[-np.pi, np.pi], [-np.pi, np.pi], [-np.pi, np.pi]],
        'distribution': ['Uniform', 'Uniform', 'Uniform']
    }
    dd_simplySupportedBeam = {
        "random variables": ['b', 'h', 'L', 'E', 'p'],
        "objective variables": ['Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8', 'Y9'],
        # 'bounds': [[0.15, 0.0075], [0.3, 0.015], [5, 0.05], [3e10, 4.5e9], [1e4, 2e3]],
        'bounds': [[0.1316, 0.1658], [0.264, 0.3323], [4.8986, 5.098], [2.4336e10, 3.8998e10], [7.5049e3, 1.5823e4]],
        'distribution': ['Lognormal', 'Lognormal', 'Lognormal', 'Lognormal', 'Lognormal']
    }
    dd_toyProblem = {
        "random variables": ['x', 'y'],
        "objective variables": ['f1', 'f2'],
        # 'bounds': [[0.15, 0.0075], [0.3, 0.015], [5, 0.05], [3e10, 4.5e9], [1e4, 2e3]],
        'bounds': [[0, 1.5], [0, 1.5]],
        'distribution': ['Uniform', 'Uniform']
    }
    dd_dqw = {
        "random variables": ['shaft_y', 'bend_out_sec_prob_length', 'bend_out_cap_gap', 'Window_margin',
                             'Shift_from_center', 'Cap_y', 'Cap_thickness', 'Cap_Height', 'Cap_Gap',
                             'Bend_out_chamber_Length', 'BP_HOM_Penetration'],
        "objective variables": ['S_0', 'f_0', 'BW'],
        'bounds': [[172.0, 258.0], [48.0, 72.0], [4.08, 6.12], [9.60, 14.40], [-18.0, -12.0], [-18.60, -12.4],
                   [2.40, 3.60], [39.68, 59.52], [12.0, 18.0], [24.40, 36.6], [16.0, 24.0]],
        'distribution': ['Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform', 'Uniform',
                         'Uniform', 'Uniform', 'Uniform']
    }

    m = UQModel()
    m.set_input_variables(dd_dqw)
    m.set_sample_size(10000)
    m.set_model(dqw_pce_from_df)
    m.run_analysis()
    # Sj, STi = m.sobol_df()
    # Sj = m.sobol_satelli()
    Sj, STj = m.sobol_janon()
    m.plot_sobol(Sj, 'Main indices', plot_type="f")
    m.plot_sobol(STj, 'Total indices')

    print(Sj)
    print()
    print(STj)
