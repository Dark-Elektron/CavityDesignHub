import math

import numpy as np
import pandas as pd
from matplotlib import text as mtext
from icecream import ic
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.stats import qmc
import itertools


class CurvedText(mtext.Text):
    """
    A text object that follows an arbitrary curve.
    """

    def __init__(self, x, y, text, axes, **kwargs):
        super(CurvedText, self).__init__(x[0], y[0], ' ', **kwargs)

        axes.add_artist(self)

        ##saving the curve:
        self.__x = x
        self.__y = y
        self.__zorder = self.get_zorder()

        ##creating the text objects
        self.__Characters = []
        for c in text:
            if c == ' ':
                ##make this an invisible 'a':
                t = mtext.Text(0, 0, 'a')
                t.set_alpha(0.0)
            else:
                t = mtext.Text(0, 0, c, **kwargs)

            # resetting unnecessary arguments
            t.set_ha('center')
            t.set_rotation(0)
            t.set_zorder(self.__zorder + 1)

            self.__Characters.append((c, t))
            axes.add_artist(t)

    # overloading some member functions, to assure correct functionality on update
    def set_zorder(self, zorder):
        super(CurvedText, self).set_zorder(zorder)
        self.__zorder = self.get_zorder()
        for c, t in self.__Characters:
            t.set_zorder(self.__zorder + 1)

    def draw(self, renderer, *args, **kwargs):
        """
        Overload of the Text.draw() function. Do not do
        do any drawing, but update the positions and rotation
        angles of self.__Characters.
        """
        self.update_positions(renderer)

    def update_positions(self, renderer):
        """
        Update positions and rotations of the individual text elements.
        """

        # preparations

        # determining the aspect ratio:
        # from https://stackoverflow.com/a/42014041/2454357

        # data limits
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()
        # Axis size on figure
        figW, figH = self.axes.get_figure().get_size_inches()
        # Ratio of display units
        _, _, w, h = self.axes.get_position().bounds
        # final aspect ratio
        aspect = ((figW * w) / (figH * h)) * (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])

        # points of the curve in figure coordinates:
        x_fig, y_fig = (
            np.array(l) for l in zip(*self.axes.transData.transform([
            (i, j) for i, j in zip(self.__x, self.__y)
        ]))
        )

        # point distances in figure coordinates
        x_fig_dist = (x_fig[1:] - x_fig[:-1])
        y_fig_dist = (y_fig[1:] - y_fig[:-1])
        r_fig_dist = np.sqrt(x_fig_dist ** 2 + y_fig_dist ** 2)

        # arc length in figure coordinates
        l_fig = np.insert(np.cumsum(r_fig_dist), 0, 0)

        # angles in figure coordinates
        rads = np.arctan2((y_fig[1:] - y_fig[:-1]), (x_fig[1:] - x_fig[:-1]))
        degs = np.rad2deg(rads)

        rel_pos = 10
        for c, t in self.__Characters:
            # finding the width of c:
            t.set_rotation(0)
            t.set_va('center')
            bbox1 = t.get_window_extent(renderer=renderer)
            w = bbox1.width
            h = bbox1.height

            # ignore all letters that don't fit:
            if rel_pos + w / 2 > l_fig[-1]:
                t.set_alpha(0.0)
                rel_pos += w
                continue

            elif c != ' ':
                t.set_alpha(1.0)

            # finding the two data points between which the horizontal
            # center point of the character will be situated
            # left and right indices:
            il = np.where(rel_pos + w / 2 >= l_fig)[0][-1]
            ir = np.where(rel_pos + w / 2 <= l_fig)[0][0]

            # if we exactly hit a data point:
            if ir == il:
                ir += 1

            # how much of the letter width was needed to find il:
            used = l_fig[il] - rel_pos
            rel_pos = l_fig[il]

            # relative distance between il and ir where the center
            # of the character will be
            fraction = (w / 2 - used) / r_fig_dist[il]

            # setting the character position in data coordinates:
            # interpolate between the two points:
            x = self.__x[il] + fraction * (self.__x[ir] - self.__x[il])
            y = self.__y[il] + fraction * (self.__y[ir] - self.__y[il])

            # getting the offset when setting correct vertical alignment
            # in data coordinates
            t.set_va(self.get_va())
            bbox2 = t.get_window_extent(renderer=renderer)

            bbox1d = self.axes.transData.inverted().transform(bbox1)
            bbox2d = self.axes.transData.inverted().transform(bbox2)
            dr = np.array(bbox2d[0] - bbox1d[0])

            # the rotation/stretch matrix
            rad = rads[il]
            rot_mat = np.array([
                [math.cos(rad), math.sin(rad) * aspect],
                [-math.sin(rad) / aspect, math.cos(rad)]
            ])

            ##computing the offset vector of the rotated character
            drp = np.dot(dr, rot_mat)

            # setting final position and rotation:
            t.set_position(np.array([x, y]) + drp)
            t.set_rotation(degs[il])

            t.set_va('center')
            t.set_ha('center')

            # updating rel_pos to right edge of character
            rel_pos += w - used


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


def load_data(filename):
    df = pd.read_excel(filename, 'Sheet1')
    return df


def plot_boxplot_multiple(df, columns_list, layout=None, xtickslabel=None):
    # set size based on number of columns
    figsize = (2.5 * len(layout[0]), 4 * len(layout))
    fig, axs = plt.subplot_mosaic(layout, figsize=figsize)

    for i, columns in enumerate(columns_list):
        all_data = df[columns]

        if i in list(xtickslabel.keys()):
            labels = xtickslabel[i]
        else:
            labels = columns

        # rectangular box plot
        bplot = axs[i].boxplot(all_data,
                               vert=True,  # vertical box alignment
                               patch_artist=True,  # fill with color
                               medianprops=dict(c='k'),
                               labels=labels)  # will be used to label x-ticks

        # fill with colors
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)

        # adding horizontal grid lines
        axs[i].yaxis.grid(True)
        # axs[i].set_ylabel('Sampled values')
        axs[i].tick_params(axis='x', labelrotation=90)
    plt.tight_layout()
    plt.show()


def plot_boxplot(df, columns, norm=False):
    if norm:
        all_data = df[columns]
        all_data_normalized = (all_data - all_data.min()) / (all_data.max() - all_data.min())
        labels = columns

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
        # fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))

        # rectangular box plot
        bplot1 = ax1.boxplot(all_data,
                             vert=True,  # vertical box alignment
                             patch_artist=True,  # fill with color
                             medianprops=dict(c='k'),
                             labels=labels)  # will be used to label x-ticks
        # ax1.set_title('Rectangular box plot')

        # plot norm
        bplot2 = ax2.boxplot(all_data_normalized,
                             vert=True,  # vertical box alignment
                             patch_artist=True,  # fill with color
                             medianprops=dict(c='k'),
                             labels=labels)  # will be used to label x-ticks
        # ax2.set_title('Normalised')

        # fill with colors
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for bplot in (bplot1, bplot2):
            for patch, color in zip(bplot['boxes'], colors):
                patch.set_facecolor(color)

        # adding horizontal grid lines
        for ax in [ax1, ax2]:
            ax.yaxis.grid(True)
            ax.set_xlabel('Geometric variables')
        ax1.set_ylabel('Sampled values')
        ax2.set_ylabel('Min-max norm. sampled values')
    else:
        all_data = df[columns]
        all_data_normalized = (all_data - all_data.min()) / (all_data.max() - all_data.min())
        labels = columns

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
        # fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))

        # rectangular box plot
        bplot = ax.boxplot(all_data,
                           vert=True,  # vertical box alignment
                           patch_artist=True,  # fill with color
                           medianprops=dict(c='k'),
                           labels=labels)  # will be used to label x-ticks

        # fill with colors
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)

        # adding horizontal grid lines
        ax.yaxis.grid(True)
        ax.set_xlabel('Geometric variables')
        ax.set_ylabel('Sampled values')

    plt.tight_layout()
    plt.show()


def calculate_discrepancy(df, lb, ub):
    space = qmc.scale(df, lb, ub, reverse=True)
    ic(qmc.discrepancy(space))


def describe_data(df):
    # print description to latex
    df_describe = df.describe()
    df_describe = df_describe.style.format(precision=2)
    print(df_describe.to_latex(position='!htb'))


def plot_spider(filename, op):
    df = pd.read_excel(filename, op).set_index('Parameter')
    df_normalised = df.div(df['2019 (CDR)'], axis=0)

    # select only particlular rows

    labels = ['Circumference [km]', 'Beam current/beam [mA]', 'Total RF voltage [GV]', 'SR loss/turn [MeV]',
              'No. of bunches/beam', 'Bunch intensity [$10^{11}$]', 'Bunch SR length [mm]', 'Bunch BS length [mm]',
              'Luminosity [$10^{34}\mathrm{cm^{-2}s^{-1}}$]']

    df_normalised = df_normalised.loc[labels]

    # Extract index values and column names for the spider chart
    categories = list(df_normalised.index)
    num_columns = len(df_normalised.columns)
    values = df_normalised.values.tolist()

    # Add the first value at the end of the list to close the chart
    values += [values[0]]
    values = np.array(values).T.tolist()

    # Create an array of evenly spaced angles for the radar chart
    angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
    angles += angles[:1]

    # Create the spider chart
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes(111, polar=True)
    labels = [l.split('[')[0] for l in labels]
    ax.set_xticks(angles[:-1], labels, size=12)
    marker = itertools.cycle(('o', 's', 'd', 'v', '^', 'p'))
    # Plot the values
    for i in range(num_columns):
        ax.plot(angles, values[i], lw=2, marker=next(marker), mec='k')

        # #adjusting plot limits
        # stretch = 0.2
        # xlim = ax.get_xlim()
        # w = xlim[1] - xlim[0]
        # ax.set_xlim([xlim[0] - stretch * w, xlim[1] + stretch * w])
        # ylim = ax.get_ylim()
        # h = ylim[1] - ylim[0]
        # ax.set_ylim([ylim[0] - stretch * h, ylim[1] + stretch * h])
        #
        # # adding the text
        # text = CurvedText(
        #     x=angles,
        #     y=values[i],
        #     text=labels[i],  # 'this this is a very, very long text',
        #     va='bottom',
        #     axes=ax,  ##calls ax.add_artist in __init__
        # )

    # ax.legend(labels=df_normalised.columns, bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left", mode="expand", ncols=6)
    ax.legend(labels=df_normalised.columns, loc="lower left", bbox_to_anchor=(1, 0))
    # Show the plot
    plt.tight_layout()
    plt.show()


def curveText(text, height, minTheta, maxTheta, ax):
    interval = np.arange(minTheta, maxTheta, .022)
    if (maxTheta <= np.pi):
        progression = interval[::-1]
        rotation = interval[::-1] - np.arctan(np.tan(np.pi / 2))
    else:
        progression = interval
        rotation = interval - np.arctan(np.tan(np.pi / 2)) - np.pi

    # Render each letter individually
    for i, rot, t in zip(progression, rotation, text):
        ax.text(i, height, t, fontsize=11, rotation=np.degrees(rot), ha='center', va='center')


if __name__ == '__main__':
    print("S T A T I S T I C S")
    # plot settings
    plot_settings()

    # pd.set_option('display.max_colwidth', None)
    # pd.set_option('display.max_columns', None)

    # filename = fr'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Data\2023_17_07_combined_results\all_results_abci_slans (including deleted results).xlsx'
    # filename = fr'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Data\all_results_abci_slans.xlsx'
    # # filename = fr'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Data\2023_17_07_combined_results\short_testing.xlsx'
    # df = load_data(filename)
    # ic(len(df.index))
    #
    # # clean dataset
    # df['freq'] = df['freq'].fillna(400.79)
    # df = df.loc[(df['freq'] >= 400.29) & (df['freq'] <= 401.29) & (df['L'] <= 200) | (df['freq'].isnull())]
    #
    # # rename columns to make shorter
    # df = df.rename(columns={'Z_long[max(0.44<f<0.77)]': 'ZL1', 'Z_long[max(0.77<f<1.1)]': 'ZL2',
    #                         'Z_long[max(1.1<f<2.0)]': 'ZL3', 'Z_trans[max(0.54<f<0.59)]': 'ZT1',
    #                         'Z_trans[max(0.59<f<0.75)]': 'ZT2', 'Z_trans[max(0.75<f<1.05)]': 'ZT3',
    #                         'Z_trans[max(1.05<f<2.0)]': 'ZT4'})
    #
    # # add new column for maximum long and transverse impedance
    # df['max_ZL'] = df[['ZL1', 'ZL2', 'ZL3']].max(axis=1)
    # df['max_ZT'] = df[['ZT1', 'ZT2', 'ZT3', 'ZT4']].max(axis=1)
    #
    # # describe dataset
    # describe_data(df[['A', 'B', 'a', 'b', 'Ri', 'L', 'Req']])
    # describe_data(df[['Epk/Eacc', 'Bpk/Eacc',
    #                   'k_loss_long', 'k_loss_trans',
    #                   'ZL1', 'ZL2', 'ZL3',
    #                   'ZT1', 'ZT2', 'ZT3', 'ZT4']])
    #
    # ic(len(df.index))
    # df = df[['A', 'B', 'a', 'b', 'Ri', 'L', 'Req', 'freq',
    #          'Epk/Eacc', 'Bpk/Eacc',
    #          'k_loss_long', 'k_loss_trans',
    #          'max_ZL', 'max_ZT',
    #          'ZL1', 'ZL2', 'ZL3',
    #          'ZT1', 'ZT2', 'ZT3', 'ZT4']]
    #
    # # plot_boxplot(df, ['A', 'B', 'a', 'b', 'Ri', 'L'])
    # # plot_boxplot(df, ['Epk/Eacc', 'Bpk/Eacc'])
    # # plot_boxplot(df, ['k_loss_long', 'k_loss_trans'])
    # # plot_boxplot(df, ['max_ZL', 'max_ZT'])
    # # plot_boxplot(df, ['ZL1', 'ZL2', 'ZL3',
    # #                   'ZT1', 'ZT2', 'ZT3', 'ZT4'])
    #
    # plot_boxplot_multiple(df,
    #                       [['A', 'B', 'a', 'b', 'Ri', 'L'], ['Epk/Eacc', 'Bpk/Eacc'], ['k_loss_long', 'k_loss_trans'],
    #                        ['ZL1', 'ZL2', 'ZL3'], ['ZT1', 'ZT2', 'ZT3', 'ZT4']],
    #                       [[0, 1, 2, 3, 4]],
    #                       xtickslabel={
    #                           0: ['$A$ [mm]', '$B$ [mm]', '$a$ [mm]', '$b$ [mm]', '$R_\mathrm{i}$ [mm]', '$L$ [mm]'],
    #                           1: ['$E_\mathrm{pk}/E_\mathrm{acc}$\n[$\cdot$]',
    #                               '$B_\mathrm{pk}/E_\mathrm{acc}$\n[mT/MV/m]'],
    #                           2: ['$k_\parallel$\n[V/pC]', '$k_\perp$\n[V/pC/m]'],
    #                           3: ['$Z_{\parallel, 1}$\n[k$\Omega$]', '$Z_{\parallel, 2}$\n[k$\Omega$]',
    #                               '$Z_{\parallel, 3}$\n[k$\Omega$]'],
    #                           4: ['$Z_{\perp, 1}$\n[k$\Omega$/m]', '$Z_{\perp, 2}$\n[k$\Omega$/m]',
    #                               '$Z_{\perp, 3}$\n[k$\Omega$/m]', '$Z_{\perp, 4}$\n[k$\Omega$/m]']})

    ############### KWT ##################

    # ################### cavity space statistics ##########
    # filename = fr'D:\Dropbox\CavityDesignHub\KWT_simulations\PostprocessingData\Data\grid_results.xlsx'
    # df = load_data(filename)
    #
    # # clean dataset
    # df['freq'] = df['freq'].fillna(801.58)
    # df = df.loc[(df['freq'] >= 801) & (df['freq'] <= 802)]
    #
    # # describe dataset
    # describe_data(df[['A_i', 'B_i', 'a_i2', 'b_i3', 'Ri_i', 'L_i', 'Req']])
    # describe_data(df[['Epk/Eacc', 'Bpk/Eacc', 'R/Q']])
    #
    # ic(len(df.index))
    # df = df[['A_i', 'B_i', 'a_i2', 'b_i3', 'Ri_i', 'L_i', 'Req', 'freq',
    #          'Epk/Eacc', 'Bpk/Eacc', 'R/Q']]
    #
    # plot_boxplot_multiple(df,
    #                       [['A_i', 'B_i', 'a_i2', 'b_i3', 'Ri_i', 'L_i'], ['Epk/Eacc', 'Bpk/Eacc', 'R/Q']],
    #                       [[0, 1]],
    #                       xtickslabel={
    #                           0: ['$A$ [mm]', '$B$ [mm]', '$a$ [mm]', '$b$ [mm]', '$R_\mathrm{i}$ [mm]', '$L$ [mm]'],
    #                           1: ['$E_\mathrm{pk}/E_\mathrm{acc}$\n[$\cdot$]',
    #                               '$B_\mathrm{pk}/E_\mathrm{acc}$\n[mT/MV/m]', '$\Omega$']})

    # ################# fcc-ee parameter evolution spider plot ###############
    # filename = fr'D:\Dropbox\CavityDesignHub\PhD_Thesis\OperatingPoints\machine_parameter_evolution.xlsx'
    # plot_spider(filename, 'W')

    ################### show variation of r/q of the first mode with cell length ####################
    # get dataframes
    filename = fr'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\modules\data_module\COMPLETE_SECOND_BATCH_9064_6D_space_w_Rsh_Q_1.xlsx'
    df1 = pd.read_excel(filename, 'Sheet1')
    filename = fr'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\modules\data_module\COMPLETE_THIRD_BATCH_9276_w_Rsh_Q_1.xlsx'
    df2 = pd.read_excel(filename, 'Sheet1')

    # get only R/Q1 and length
    df1_ = df1.loc[:, ['L', 'Rsh/Q1', 'Rsh/Q']]
    df2_ = df2.loc[:, ['L', 'Rsh/Q1', 'Rsh/Q']]

    # combine dataframe
    data = pd.concat([df1_, df2_], axis=0)
    ic(len(data.index))

    # Define the bin ranges for variable 'x'
    bins = [140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190]

    # Use cut to categorize 'x' into bins
    labels = [f'{bins[i]}-{bins[i + 1]}' for i in range(len(bins) - 1)]
    data['L_bins'] = pd.cut(data['L'], bins=bins, labels=labels)

    # Calculate the average of 'y' in each bin
    bin_averages = data.groupby('L_bins')['Rsh/Q1', 'Rsh/Q'].mean()

    # Create a bar plot
    # bin_averages.plot(kind='bar')
    # Create a box plot

    fig, ax = plt.subplot_mosaic([[0], [1]], sharex=True, figsize=(12, 4))
    ax_ylabels = ['$R/Q_\mathrm{TM010-0}$ [$\Omega$]', '$R/Q_\mathrm{TM010-\pi}$ [$\Omega$]']
    for i, pa in enumerate(['Rsh/Q1', 'Rsh/Q']):
        ax[i].plot(data['L'], data[pa], ls='', marker='o', ms=3, c='gray', mfc='none', alpha=0.1)
        ax[i].boxplot([data[data['L_bins'] == bin][pa] for bin in labels], labels=labels,
                    vert=True,  # vertical box alignment
                    patch_artist=True,  # fill with color
                    medianprops=dict(c='k'),
                    positions=[(bins[i]+bins[i+1])/2 for i in range(len(bins) - 1)])
        ax[i].set_ylabel(ax_ylabels[i])
    # Label the x-axis with bin ranges
    plt.xlabel('L [mm]')

    # Label the y-axis as needed
    # plt.ylabel('$R/Q$ [$\Omega$]')
    plt.tight_layout()
    # Show the plot
    plt.show()
