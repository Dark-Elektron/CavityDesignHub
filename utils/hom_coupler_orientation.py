import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate
from utils.file_reader import FileReader
import matplotlib as mpl
import seaborn as sns
fr = FileReader()

def plot_settings():
    mpl.rcParams['xtick.labelsize'] = 18
    mpl.rcParams['ytick.labelsize'] = 18

    mpl.rcParams['axes.labelsize'] = 20
    mpl.rcParams['axes.titlesize'] = 20
    mpl.rcParams['legend.fontsize'] = 20
    mpl.rcParams['legend.title_fontsize'] = 20

    mpl.rcParams['figure.figsize'] = [10, 6]
    mpl.rcParams['figure.dpi'] = 100



# file = fr.excel_reader(r"D:\Dropbox\HC_FPC_Configs.xlsx")
# print(file["2hc_2fpc"]["angle_corr"].to_list())
# categories = file["2hc_2fpc"]["angle_corr"].to_list()
#
# angle = file["2hc_2fpc"]["angle_corr"].to_list()
# hc2_fpc1 = file["2hc_1fpc"]["Z"].to_list()
# hc2_fpc2 = file["2hc_2fpc"]["Z"].to_list()
#
# plt.rcParams.update({
#     "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # red   with alpha = 30%
#     "axes.facecolor":    (1.0, 1.0, 1.0, 0.0),  # green with alpha = 50%
#     "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # blue  with alpha = 20%
# })
# fig, ax = plt.subplots(figsize=(9.8, 6)) #, subplot_kw = dict(polar = True)
# ax.plot(angle, hc2_fpc1, label='2HC1FPC', marker='o', lw=3)
# ax.plot(angle, hc2_fpc2, label='2HC2FPC', marker='o', lw=3)
# ax.set_xlabel(r"$\alpha$")
# ax.set_ylabel(r"Z_$_{\mathrm{T, TM_{110}}} \mathrm{k\Omega/m}$")
# ax.legend()
# plt.tight_layout()
# plt.show()


if __name__ == '__main__':
    plot_settings()
    #####################
    from matplotlib import ticker, cm, colors
    # surface plot
    file = fr.excel_reader(r"D:\Dropbox\HC_FPC_Configs.xlsx")
    sheet = "4hc_1fpc_2angles"
    ang = (file[sheet]["angle"]).to_list()
    ang_fpc = file[sheet]["angle_fpc"].to_list()

    X = np.reshape(ang, (10, 14))
    Y = np.reshape(ang_fpc, (10, 14))
    Z = file[sheet]["ZT"].to_list()
    Z2 = np.reshape(Z, np.shape(X))

    xnew, ynew, znew = X, Y, Z2
    print(xnew)
    print(ynew)

    fig = plt.figure(figsize=(18, 4))
    ax = fig.gca()
    xticks, yticks = xnew[0], ynew.T[0]
    g = sns.heatmap(znew, cmap=cm.Pastel1, annot=True, norm=colors.LogNorm(vmin=10, vmax=25), linewidths=.5, ax=ax, xticklabels=xticks, yticklabels=yticks,
                    cbar_kws=dict(pad=0.01), annot_kws={"size": 14})
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')

    plt.xlabel(r"$\beta [\mathrm{^\circ}]$")
    plt.ylabel(r"$\alpha [\mathrm{^\circ}]$")
    fig.tight_layout()

    plt.show()
