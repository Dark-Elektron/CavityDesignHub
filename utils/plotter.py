
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
# from numpy._distributor_init import Book1


class Plotter:
    def __init__(self):
        self.color = ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c']

    def formatAxis(self, ax, color, position, offset=None):
        if offset is None:
            offset = {'right': 0}

        for pos in position:
            ax.spines[pos].set_color(color)

        for key, val in offset.items():
            ax.spines[key].set_position(('outward', val))  # change from outward to a more generalised keyword later

        ax.tick_params(axis='x', colors=color)
        ax.tick_params(axis='y', colors=color)

    def plot(self, sheet_name, x_name, y_name, ax1_name, ax2_name, ax3_name):
        file = pd.ExcelFile(r'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Dataframe_full_data.xlsx')
        print(file.sheet_names)

        dd = {}
        dd[sheet_name] = file.parse(sheet_name)
        # for sht_name in file.sheet_names:
        #     dd[sheet_name] = file.parse(sht_name)

        print(dd[sheet_name])

        # reshape data
        x_data = np.array(dd[sheet_name][x_name])[0:70].reshape(10, 7) # take first multiple of 7 elements
        y_data = np.array(dd[sheet_name][y_name][0:70]).reshape(10, 7)
        ax1_data = np.array(dd[sheet_name][ax1_name])[0:70].reshape(10, 7)
        ax2_data = np.array(dd[sheet_name][ax2_name])[0:70].reshape(10, 7)
        ax3_data = np.array(dd[sheet_name][ax3_name])[0:70].reshape(10, 7)
        # print(f'xdata: {x_data}')

        # create figure and axis objects with subplots()
        fig, ax = plt.subplots()
        ls = ['-', '--', '-.', ':']
        mk = ['o', 'D', 'p', 'X']
        for i, x in enumerate(x_data):
            # make a plot
            ax.plot(x, y_data[i], color=self.color[0], linestyle=ls[0], marker=mk[0])
            # ax.scatter(x, y_data[i], color=self.color[0], linestyle=ls[0], marker=mk[0])
            # set x-axis label
            ax.set_xlabel(x_name, fontsize=14)
            # set y-axis label
            ax.set_ylabel(y_name, color=self.color[0], fontsize=14)
            ax.tick_params(axis='y', colors=self.color[0])

            # SECOND PLOT
            # twin object for two different y-axis on the sample plot
            if i == 0:
                ax1 = ax.twinx()
                ax1.set_ylabel(ax1_name, color=self.color[1], fontsize=14)
                self.formatAxis(ax1, self.color[1], ['right'])
            # make a plot with different y-axis using second axis object
            ax1.plot(x, ax1_data[i], color=self.color[1], linestyle=ls[1], marker=mk[1])
            # ax1.scatter(x, ax1_data[i], color=self.color[1], linestyle=ls[1], marker=mk[1])

            # twin object for two different y-axis on the sample plot
            if i == 0:
                ax2 = ax.twinx()
                ax2.set_ylabel(ax2_name, color=self.color[2], fontsize=14)
                self.formatAxis(ax2, self.color[i+2], ['right'], offset={'right': 60})
            # make a plot with different y-axis using second axis object
            ax2.plot(x, ax2_data[i], color=self.color[2], linestyle=ls[2], marker=mk[2])
            # ax2.scatter(x, ax2_data[i], color=self.color[2], linestyle=ls[2], marker=mk[2])

            # twin object for two different y-axis on the sample plot
            if i == 0:
                ax3 = ax.twinx()
                ax3.set_ylabel(ax3_name, color=self.color[3], fontsize=14)
                self.formatAxis(ax3, self.color[3], ['right'], offset={'right': 2*60})
            # make a plot with different y-axis using second axis object
            ax3.plot(x, ax3_data[i], color=self.color[3], linestyle=ls[3], marker=mk[3])
            # ax3.scatter(x, ax3_data[i], color=self.color[3], linestyle=ls[3], marker=mk[3])

        fig.tight_layout()
        plt.show()

        # save the plot as a file
        # fig.savefig('two_different_y_axis_for_single_python_plot_with_twinx.jpg',
        #             format='jpeg',
        #             dpi=100,
        #             bbox_inches='tight')


if __name__ == '__main__':
    plotter = Plotter()
    plotter.plot('Midcell node_editor', 'Ri', 'ZTM011', 'ZTM012', 'ZD1', 'ZD2')
