import matplotlib
import numpy as np
import oapackage
import pandas as pd
from icecream import ic
from matplotlib import pyplot as plt, animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import matplotlib as mpl


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


def pareto_front(df, columns, ax=None, show='all', label='', kwargs_dataframe=None, kwargs_pareto=None):
    if kwargs_dataframe is None:
        kwargs_dataframe = {}
    if kwargs_pareto is None:
        kwargs_pareto = {}

    datapoints = df.loc[:, columns] * (-1)
    pareto = oapackage.ParetoDoubleLong()
    # ic(datapoints)
    for ii in range(0, datapoints.shape[0]):
        w = oapackage.doubleVector(tuple(datapoints.iloc[ii].values))
        pareto.addvalue(w, ii)
    # pareto.show(verbose=1)  # Prints out the results from pareto

    lst = pareto.allindices()  # the indices of the Pareto optimal designs
    poc = len(lst)  # number of pareto shapes
    reorder_idx = list(lst) + [i for i in range(len(df)) if i not in lst]  # reordered index putting pareto shapes first

    pareto_shapes = df.loc[lst, :]
    # sort pareto shapes in ascending x axis data
    pareto_shapes = pareto_shapes.sort_values(by=[columns[0]])

    if ax:
        if show == 'all':
            ax.scatter(df[columns[0]], df[columns[1]], **kwargs_dataframe)

        ax.plot(pareto_shapes[columns[0]], pareto_shapes[columns[1]], **kwargs_pareto, label=label)

    return pareto_shapes


def plot_pareto_surface(df, columns, ax, kwargs=None):
    if kwargs is None:
        kwargs = {}
    pareto = pareto_front(df, columns)

    x, y, z = pareto['Epk/Eacc'], pareto['Bpk/Eacc'], pareto['R/Q']
    xi, yi = np.meshgrid(np.linspace(min(x), max(x), 100),
                         np.linspace(min(y), max(y), 100))
    zi = griddata((x, y), z, (xi, yi), method='cubic')
    surf = ax.plot_surface(xi, yi, zi, antialiased=False, **kwargs)

    return surf


def create_evolution_animation(frame):
    for ax_index in axs:
        axs[ax_index].clear()

    ################### plot 2d ##################################
    grid_results_folder = fr'D:\Dropbox\CavityDesignHub\KWT_simulations\PostprocessingData\Data\grid_results.xlsx'

    columns_array = [['Epk/Eacc', 'Bpk/Eacc'], ['Epk/Eacc', 'R/Q'], ['Bpk/Eacc', 'R/Q']]
    columns_par_array = [['A', 'B'], ['a', 'b'], ['A', 'Ri']]
    cmap = matplotlib.colormaps['Pastel2_r']
    norm = matplotlib.colors.Normalize(vmin=0, vmax=49)

    for i, (columns, columns_par) in enumerate(zip(columns_array, columns_par_array)):
        grid_results = pd.read_excel(grid_results_folder, 'Sheet1')
        qmc_pareto_shapes = pareto_front(grid_results, columns, axs[i], show='pareto',
                                         label=f'QMC \n({len(grid_results.index)} geoems.)',
                                         kwargs_dataframe={'facecolors': 'none', 'edgecolor': 'grey'},
                                         kwargs_pareto={'c': 'k', 'marker': 'o', 'mec': 'k'})

        for r in [frame]:
            for lab, opt_code_result_folder, error_file in zip(['LHS2'],
                                                               [
                                                                   fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS\Generation{r}.xlsx'],
                                                               [
                                                                   fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS_LHS2\inerp_error_and_average.txt']):
                opt_code_result = pd.read_excel(opt_code_result_folder, 'Sheet1')

                pareto_shapes = pareto_front(opt_code_result, columns, axs[i], show='pareto',
                                             label=f"{lab}: G{r} \n({len(opt_code_result.index)} geoms).",
                                             kwargs_dataframe={'facecolors': 'none', 'edgecolor': 'b'},
                                             kwargs_pareto={'marker': 'o',
                                                            'mec': 'k'},
                                             # kwargs_pareto={'c': f'{matplotlib.colors.rgb2hex(cmap(norm(r)))}', 'marker': 'o',
                                             #                'mec': 'k'}
                                             )
                axs[i + 3].scatter(grid_results[columns_par[0]], grid_results[columns_par[1]], s=5)
                axs[i + 3].scatter(opt_code_result[columns_par[0]], opt_code_result[columns_par[1]], s=5)
                axs[i + 3].scatter(qmc_pareto_shapes[columns_par[0]], qmc_pareto_shapes[columns_par[1]], c='r',
                                   label='qmc', s=5)
                axs[i + 3].scatter(pareto_shapes[columns_par[0]], pareto_shapes[columns_par[1]], c='b', edgecolor='k',
                                   label='ea', s=5)

                # load interpolation error
                error = pd.read_csv(error_file, header=None, sep='\s+')
                # axs[3].plot(error[0], label=lab)
                # axs[4].plot(error[1], label=lab)

        axs[i].set_xlabel(columns[0])
        axs[i].set_ylabel(columns[1])
        axs[i + 3].set_xlabel(columns_par[0])
        axs[i + 3].set_ylabel(columns_par[1])

    # axs[i].legend(bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=10, loc='lower left', mode='expand')
    axs[0].legend()
    axs[3].legend()

    # axs[3].set_yscale('log')
    # axs[4].set_yscale('log')
    # axs[3].set_xlabel('Interpolation error')
    # axs[3].set_ylabel()
    # plot error

    # plt.tight_layout()
    # plt.show()


if __name__ == '__main__':
    plot_settings()

    ################## plot 2d ##################################
    grid_results_folder = fr'D:\Dropbox\CavityDesignHub\KWT_simulations\PostprocessingData\Data\grid_results.xlsx'
    fig, axs = plt.subplot_mosaic([[0, 1, 2], [3, 4, 5]], figsize=(12, 6))

    columns_array = [['Epk/Eacc', 'Bpk/Eacc'], ['Epk/Eacc', 'R/Q'], ['Bpk/Eacc', 'R/Q']]
    columns_array = [['Epk/Eacc []', 'Bpk/Eacc [mT/MV/m]']]
    columns_par_array = [['A', 'B'], ['a', 'b'], ['A', 'Ri']]
    cmap = matplotlib.colormaps['Pastel2_r']
    norm = matplotlib.colors.Normalize(vmin=0, vmax=49)
    ff = fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData'
    for i, (columns, columns_par) in enumerate(zip(columns_array, columns_par_array)):
        # grid_results = pd.read_excel(grid_results_folder, 'Sheet1')
        # qmc_pareto_shapes = pareto_front(grid_results, columns, axs[i], show='pareto',
        #                                  label=f'MC (n={len(grid_results.index)})',
        #                                  kwargs_dataframe={'facecolors': 'none', 'edgecolor': 'grey'},
        #                                  kwargs_pareto={'c': 'k', 'marker': 'o', 'mec': 'k', 'ms': 3})

        for r in range(20):
            # for lab, opt_code_result_folder, error_file in zip(['LHS', 'LHS2', 'Random'],
            #                                                    [fr'{ff}\SLANS_LHS\Generation{r}.xlsx',
            #                                                     fr'{ff}\SLANS_LHS2\Generation{r}.xlsx',
            #                                                     fr'{ff}\SLANS_kwt_random1\Generation{r}.xlsx'],
            #                                                    [fr'{ff}\SLANS_LHS2\inerp_error_and_average.txt',
            #                                                     fr'{ff}\SLANS_LHS2\inerp_error_and_average.txt',
            #                                                     fr'{ff}\SLANS_LHS2\inerp_error_and_average.txt']):

            for lab, opt_code_result_folder in zip(['LHS'], [fr'D:\Dropbox\CavityDesignHub\MuCol_Study\SimulationData\ConsoleTest\cavities\SimulationData\Optimisation\Generation{r}.xlsx']):
                opt_code_result = pd.read_excel(opt_code_result_folder, 'Sheet1')

                pareto_shapes = pareto_front(opt_code_result, columns, axs[i], show='none',
                                             label=f"{lab}: g{r} ($n$={len(opt_code_result.index)}).",
                                             kwargs_dataframe={'facecolors': 'none', 'edgecolor': 'b'},
                                             kwargs_pareto={'marker': 'o', 'mec': 'k', 'ms': 3},
                                             # kwargs_pareto={'c': f'{matplotlib.colors.rgb2hex(cmap(norm(r)))}', 'marker': 'o',
                                             #                'mec': 'k'}
                                             )
                # axs[i+3].plot(grid_results[columns[0]], grid_results[columns[1]], marker='o', ms=5, lw=0)
                axs[i+3].plot(opt_code_result[columns[0]], opt_code_result[columns[1]], marker='o', ms=5, lw=0)
                # axs[i+3].plot(qmc_pareto_shapes[columns[0]], qmc_pareto_shapes[columns[1]], marker='o', c='r', label='qmc', ms=5, lw=0)
                axs[i+3].plot(pareto_shapes[columns[0]], pareto_shapes[columns[1]], marker='o', c='b', mec='k', label='ea', ms=5, lw=0)

                # # load interpolation error
                # error = pd.read_csv(error_file, header=None, sep='\s+')
                # axs[3].plot(error[0], label=lab)
                # axs[4].plot(error[1], label=lab)


        axs[i].set_xlabel(columns[0])
        axs[i].set_ylabel(columns[1])
        axs[i+3].set_xlabel(columns_par[0])
        axs[i+3].set_ylabel(columns_par[1])
    lines, labels = axs[0].get_legend_handles_labels()
    # axs[i].legend(bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=10, loc='lower left', mode='expand')
    fig.legend(*axs[0].get_legend_handles_labels(), loc="upper left", mode="expand", ncol=4)
    axs[3].legend()

    axs[3].set_yscale('log')
    axs[4].set_yscale('log')
    axs[3].set_xlabel('Interpolation error')
    # axs[3].set_ylabel()
    # plot error

    plt.tight_layout()
    plt.show()

    # ################### plot surface #########################
    # grid_results_folder = fr'D:\Dropbox\CavityDesignHub\KWT_simulations\PostprocessingData\Data\grid_results.xlsx'
    # opt_code_result_folder_lhs = fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS_opt_KWT\Generation49.xlsx'
    # opt_code_result_folder_random = fr'D:\Dropbox\CavityDesignHub\Cavity800\SimulationData\SLANS\Generation49.xlsx'
    #
    # grid_results = pd.read_excel(grid_results_folder, 'Sheet1')
    # opt_code_result_lhs = pd.read_excel(opt_code_result_folder_lhs, 'Sheet1')
    # opt_code_result_random = pd.read_excel(opt_code_result_folder_random, 'Sheet1')
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    # ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    # # ax3 = fig.add_subplot(1, 2, 3, projection='3d')
    # axs = [ax1, ax2]
    #
    # grid_pareto_surface = plot_pareto_surface(grid_results, ['Epk/Eacc', 'Bpk/Eacc', 'R/Q'], axs[0],
    #                                           {'cmap': 'gray', 'edgecolor': 'k'})
    # opt_pareto_surface_lhs = plot_pareto_surface(opt_code_result_lhs, ['Epk/Eacc', 'Bpk/Eacc', 'R/Q'], axs[1],
    #                                              {'cmap': 'RdBu', 'edgecolor': 'k'})
    # # opt_pareto_surface_random = plot_pareto_surface(opt_code_result_random, ['Epk/Eacc', 'Bpk/Eacc', 'R/Q'], axs[2],
    # #                                                 {'cmap': 'RdBu', 'edgecolor': 'k'})
    #
    # plt.show()

    # fig, axs = plt.subplot_mosaic([[0, 1, 2], [3, 4, 5]], figsize=(12, 5.5))
    # ani = animation.FuncAnimation(fig=fig, func=create_evolution_animation, frames=49, interval=1000)
    # ani.save(filename="D:\Dropbox\Quick presentation files/ffmpeg_example.gif", writer="pillow")
    # # plt.show()
