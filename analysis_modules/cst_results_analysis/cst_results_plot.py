import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from icecream import ic
from matplotlib import image as mpimg
from scipy.interpolate import make_interp_spline, BSpline


def plot_settings():
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12

    mpl.rcParams['axes.labelsize'] = 12
    mpl.rcParams['axes.titlesize'] = 12
    mpl.rcParams['legend.fontsize'] = 12
    mpl.rcParams['legend.title_fontsize'] = 12

    mpl.rcParams['figure.dpi'] = 100

    # Set the desired colormap
    plt.rcParams['axes.prop_cycle'] = plt.cycler('color', plt.cm.Set2.colors)


def interpolate_spline(x, y, n=500):
    # 300 represents number of points to make between T.min and T.max
    xnew = np.linspace(x.min(), x.max(), n)

    spl = make_interp_spline(x, y, k=3)
    y_smooth = spl(xnew)

    return xnew, y_smooth


def generic_1d_plot(file_folder, layout=None, sep='\s+', labels=None, other_info=None, fig=None, header=None,
                    logx=False, logy=False, legend=True, smooth=False, line_kwargs=None, yscale='linear'):
    if line_kwargs is None:
        line_kwargs = {}
    axs_already_defined = False
    if labels is None:
        labels = {}
    if other_info is None:
        other_info = {}

    if fig is None:
        figsize = (4.5 * len(layout[0]), 3.5 * len(layout))
        fig, axs = plt.subplot_mosaic(layout, figsize=figsize)
    else:
        axs = fig.axes
        axs_already_defined = True

    for n, file in enumerate(file_folder):
        if n in list(other_info.keys()):
            if 'kind' in list(other_info[n].keys()):
                if other_info[0]['kind'] == 'image':
                    img = mpimg.imread(file)
                    axs[n].imshow(img)
                    axs[n].axis('off')
        else:
            if header is None:
                df = pd.read_csv(file, sep=sep, header=None)
            else:
                df = pd.read_csv(file, sep=sep)

            default_labels = list(df.columns)
            x, y = df[default_labels[0]], df[default_labels[1]]
            if yscale == 'dB':
                y = 10*np.log10(np.sqrt(df[default_labels[1]]**2 + df[default_labels[2]]**2))

            if logy and not logx:
                # interpolate spline
                if smooth:
                    x, y = interpolate_spline(x, np.log(y), 10000)
                axs[n].semilogy(x, y, **line_kwargs)
            elif logx and not logy:
                # interpolate spline
                if smooth:
                    x, y = interpolate_spline(np.log(x), y, 10000)
                axs[n].semilogx(x, y, **line_kwargs)
            elif not logx and not logy:
                # interpolate spline
                if smooth:
                    x, y = interpolate_spline(x, y, 10000)
                axs[n].plot(x, y, **line_kwargs)
            else:
                # interpolate spline
                if smooth:
                    x, y = interpolate_spline(np.log(x), np.log(y), 10000)
                axs[n].loglog(x, y, **line_kwargs)

            axs[n].set_xlim(min(df[default_labels[0]]))

            # if n in list(columns.keys()):
            #     axs[n].set_xlabel(columns[n][0])
            #     axs[n].set_ylabel(columns[n][1])

            if n in list(labels.keys()):
                axs[n].set_xlabel(labels[n][0])
                axs[n].set_ylabel(labels[n][1])
                axs[n].lines[-1].set_label(labels[n][2])
            else:
                if not axs_already_defined:
                    axs[n].set_xlabel(default_labels[0])
                    axs[n].set_ylabel(default_labels[1])

    return fig


def plot_field_value_along_curve(file_folder, cst_filenames, results, layout, labels):
    fig, axs = plt.subplot_mosaic(layout, figsize=(8, 4))
    for n, cst_filename in enumerate(cst_filenames):
        filename = [fr'{file_folder}\{cst_filename}\Export\{x}.txt' for x in results]

        generic_1d_plot(filename, layout, labels=labels[n], fig=fig, logy=True)
    lines_labels = axs[0].get_legend_handles_labels()
    lines, labels = lines_labels
    fig.axes[0].legend(bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=3, loc='lower left', mode='expand')
    # fig.legend(lines, labels, loc="outside upper center", mode="expand", ncol=6)
    # fig.subplots_adjust(hspace=0.4)
    plt.tight_layout()
    return fig


def kb_vs_lfd(file_folder):
    file_folder = [fr'{file_folder}\k_lfd_vs_kb.txt']
    fig = generic_1d_plot(file_folder, [[0]],
                          labels={0: ['$k_\mathrm{b}$ [MN/m]', '$K_\mathrm{LFD}~\mathrm{[Hz/(MV/m)^2]}$',
                                      r'$K_\mathrm{LFD}$']
                                  }, header=True, logx=True, logy=False, legend=False, smooth=False,
                          line_kwargs={'marker': 'o', 'mec': 'k'})
    fig.axes[0].legend()
    plt.tight_layout()
    return fig


def kb_vs_pressure(file_folder):
    file_folder = [fr'{file_folder}\k_pressure_vs_kb.txt']
    fig = generic_1d_plot(file_folder, [[0]],
                          labels={0: ['$k_\mathrm{b}$ [MN/m]', '$K_\mathrm{p}$ [Hz/mbar]', r'$K_\mathrm{p}$']
                                  }, header=True, logx=True, logy=False, legend=False, smooth=False,
                          line_kwargs={'marker': 'o', 'mec': 'k'})
    fig.axes[0].legend()
    plt.tight_layout()
    return fig


if __name__ == '__main__':
    plot_settings()
    # w_op_folder, ttbar_op_folder = fr'D:\CST Studio\3. W\Mechanical', fr'D:\CST Studio\5. tt\Mechanical'
    #
    # kb = [0.4509, 1.12725, 2.2545, 3.38175, 4.509, 22.545, 45.09, 450.9, 'fe']
    # cst_files_lfd = [fr'SM_C3794_kb_{x}_LFD' if isinstance(x, float) else fr'SM_C3794_fe_LFD' for x in kb]
    # cst_files_press = [fr'SM_C3794_kb_{x}_Pressure' if isinstance(x, float) else fr'SM_C3794_fe_Pressure' for x in kb]
    #
    # layout = [[0]]
    # res = ['curve1_Displacement']
    # labels = [{0: ['$L$ [mm]', 'Displacement [mm]', '$k_\mathrm{b} =' + fr'{x}' + '~\mathrm{MN/m}$']}
    #           if isinstance(x, float) else {0: ['$L$ [mm]', 'Displacement [mm]', 'Fixed ends']} for x in kb]

    # res = ['curve1_Von Mises']
    # labels = [{0: ['$L$ [mm]', 'Von Mises [GPa]', '$k_\mathrm{b} ='+ fr'{x}' + '~\mathrm{MN/m}$']}
    #               if isinstance(x, float) else {0: ['$L$ [mm]', 'Von Mises [GPa]', 'Fixed ends']} for x in kb]

    # w
    # fig = plot_field_value_along_curve(w_op_folder, cst_files_lfd, results=res, layout=layout, labels=labels)
    # fig = plot_field_value_along_curve(w_op_folder, cst_files_press, results=res, layout=layout, labels=labels)
    # fig = kb_vs_lfd(w_op_folder)
    # fig = kb_vs_pressure(w_op_folder)

    # # ttbar
    # cst_files_lfd = [fr'SM_C3795_kb_{x}_LFD' if isinstance(x, float) else fr'SM_C3795_fixed_ends_LFD' for x in kb]
    # cst_files_press = [fr'SM_C3795_kb_{x}_Pressure' if isinstance(x, float) else fr'SM_C3795_fixed_ends' for x in kb]
    #
    # fig = plot_field_value_along_curve(ttbar_op_folder, cst_files_lfd, results=res, layout=layout, labels=labels)
    # fig = plot_field_value_along_curve(ttbar_op_folder, cst_files_press, results=res, layout=layout, labels=labels)
    # # fig = kb_vs_lfd(ttbar_op_folder)
    # # fig = kb_vs_pressure(ttbar_op_folder)
    #
    # plt.show()
    #
    # df_in = pd.read_csv(fr'D:\CST Studio\3. W\Wakefield\W_C3794\Export\Power_Excitation (pb)_Power Accepted.txt',
    #                     sep='\s+', header=None)
    # df_port_1 = pd.read_csv(
    #     fr'D:\CST Studio\3. W\Wakefield\W_C3794\Export\Power_Excitation (pb)_Power Accepted per Port_Port 1.txt',
    #     sep='\s+', header=None)
    # df_port_2 = pd.read_csv(
    #     fr'D:\CST Studio\3. W\Wakefield\W_C3794\Export\Power_Excitation (pb)_Power Accepted per Port_Port 2.txt',
    #     sep='\s+', header=None)
    # ic(np.sum(df_in[1]))
    # ic(np.sum(df_port_1[1]) / np.sum(df_in[1]), np.sum(df_port_2[1]) / np.sum(df_in[1]))

    fig, axs = plt.subplot_mosaic([[0]], figsize=(10, 3))
    filefolders = [r'D:\CST Studio\3. W\Couplers\QRWG3\Export\RWG.txt', r'D:\CST Studio\3. W\Couplers\QRWG3\Export\QRWG.txt']
    labels= [{'0': ['RWG']}, {'0': ['QRWG']}]
    for i, (folder, labels) in enumerate(zip(filefolders, labels)):
        generic_1d_plot([folder], fig=fig, labels=labels, yscale='dB')

    plt.legend()
    plt.minorticks_on()
    plt.show()
