import json
import math
import os
import re
import shutil
import time
import matplotlib
import pandas as pd
from PyQt5.QtWidgets import QMessageBox

from utils.file_reader import FileReader
import matplotlib.pyplot as plt
import matplotlib as mpl
# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']  # for \text command
from labellines import labelLine, labelLines
import multiprocessing as mp
fr = FileReader()

import scipy.signal as sps
import numpy as np
from termcolor import colored

file_color = 'cyan'

# DEBUG = False
DEBUG = True


def print_(*arg):
    if DEBUG: print(colored(f'\t\t{arg}', file_color))


class ABCIData:
    def __init__(self, dir, fid, MROT):
        self.title_dict = {}
        self.data_dict = {}
        self.dir = dir
        self.fid = fid
        self.MROT = MROT
        self.loss_factor = {}
        self.SIG = None
        self.wakelength = None
        self.checks(dir, fid, MROT)

    def checks(self, dirc, fid, MROT):
        dirc = f'{dirc}/{fid}/Cavity_MROT_{MROT}.top'
        print(dirc)
        if os.path.exists(dirc):
            self._get_plot_data(dirc)
        else:
            print_("Hey chief, there seems to be a problem with the file directory. Please check.")

    def _get_plot_data(self, dir):
        frame_objects = {}
        frame_titles_objects = {}
        frame_count = 0

        plot_decorations = [r'Cavity Shape Input', r'Cavity Shape Used', r'Wake Potentials',
                            r'Real Part of Longitudinal Impedance', r'Imaginary Part of Longitudinal Impedance',
                            'Frequency Spectrum of Loss Factor', r'Loss Factor Spectrum Integrated upto F',
                            r'Real Part of Long. + Log Impedance', r'Imaginary Part of Long. + Log Impedance',
                            r'Spectrum of Long. + Log Loss Factor', r'Long. + Log Factor Integrated up to F',
                            r'Real Part of Azimuthal Impedance',
                            r'Imaginary Part of Azimuthal Impedance', r'Real Part of Transverse Impedance',
                            r'Imaginary Part of Transverse Impedance']

        with open(dir, 'r') as f:
            for line in f:
                if 'NEW FRAME' in line:
                    # add frame to frame_objects
                    if frame_count > 0:
                        frame_titles_objects[frame_count - 1] = frame_title

                        # select suitable title
                        for decor in plot_decorations:
                            if decor in frame_title[0] + frame_title[1]:
                                frame_objects[decor] = line_objects
                                break

                    # reset lists
                    line_objects = []
                    frame_title = []
                    new_frame = True

                    frame_count += 1

                if new_frame:
                    new_frame = False

                # add titles to frame_title
                if 'TITLE' in line or 'MORE' in line:
                    frame_title.append(line)
                    # print(line)
                    # if 'Transverse Wake' in line:
                    #     key = 'Transverse'
                    #     print('True')

                    # get other parameters from title
                    try:
                        if 'Azimuthal Wake' in line:
                            key = 'Azimuthal'
                        if 'Transverse Wake' in line:
                            key = 'Transverse'
                        if 'Longitudinal Wake' in line:
                            key = 'Longitudinal'

                        if 'Loss Factor' in line and key:
                            indx_0 = line.index('= ')
                            indx_1 = line.index('V/pC')
                            loss_factor = float(line[indx_0 + 2:indx_1])

                            self.loss_factor[key] = loss_factor
                            key = None

                        # get sigma
                        if 'SIG' in line and self.SIG is None:
                            indx_0 = line.index('SIG= ')
                            indx_1 = line.index(' cm')
                            self.SIG = float(line[indx_0 + 5:indx_1])
                    except:
                        pass

                if 'SET LIMITS' in line:
                    x, y = [], []

                    # get wakelength
                    if 'SET LIMITS X   0.00000E+00' in line:
                        indx_0 = line.index('SET LIMITS X   0.00000E+00')
                        indx_1 = line.index('   Y')
                        self.wakelength = float(line[indx_0 + len('SET LIMITS X   0.00000E+00'):indx_1])

                if 'JOIN' in line:
                    new_line_object = True
                else:
                    new_line_object = False

                if new_line_object:
                    # add x and y to line objects
                    line_objects.append([x, y])

                    # reset x and y
                    x, y = [], []

                if len(line.strip().split()) == 2 and line.strip().split()[0] != 'NEW' and line.strip().split()[
                    0] != 'JOIN':
                    new_line_object = False
                    ll = [float(a) for a in line.strip().split()]
                    x.append(ll[0])
                    y.append(ll[1])

            # EOF append last list
            frame_titles_objects[frame_count - 1] = frame_title
            # select suitable title
            for decor in plot_decorations:
                if decor in frame_title[0] + frame_title[1]:
                    frame_objects[decor] = line_objects
                    break

            self.data_dict = {}
            for key, lines in frame_objects.items():
                x, y = [], []
                for line in lines:
                    if len(line[0]) > 2:
                        x += line[0]
                        y += line[1]

                self.data_dict[key] = [x, y]

    def get_data(self, key):
        plot_decorations = [r'Cavity Shape Input', r'Cavity Shape Used', r'Wake Potentials',
                            r'Real Part of Longitudinal Impedance', r'Imaginary Part of Longitudinal Impedance',
                            'Frequency Spectrum of Loss Factor', r'Loss Factor Spectrum Integrated upto F',
                            r'Real Part of Long. + Log Impedance', r'Imaginary Part of Long. + Log Impedance',
                            r'Spectrum of Long. + Log Loss Factor', r'Long. + Log Factor Integrated up to F',
                            r'Real Part of Azimuthal Impedance',
                            r'Imaginary Part of Azimuthal Impedance', r'Real Part of Transverse Impedance',
                            r'Imaginary Part of Transverse Impedance']
        if isinstance(key, int):
            key = plot_decorations[key]

        # print_(self.data_dict)
        self.x = self.data_dict[key][0]
        self.y = self.data_dict[key][1]

        # Process all needed data
        self.x_peaks, self.y_peaks = self.get_peaks(threshold=max(self.y) * 0.02)

        # if key != 'Loss Factor Spectrum Integrated up to F':
        #     self.bands = self._get_bands()
        #     self.fl_list = self._interp_fl()
        #     self.fr_list = self._interp_fr()
        #     self.BW = self._get_BW()
        #     self.Q = self._get_Q()
        #     self.R_Q = self._get_RQ()

        return self.x, self.y, self.fid

    def _get_title(self, MROT):

        if MROT == 0:
            self.plot_decorations = {0: [r'Cavity Shape Input', r'Z-axis (m)', 'R-axis (m)'],
                                     1: [r'Cavity Shape Used', r'Z-axis (m)', 'R-axis (m)'],
                                     2: [r'Wake Potentials', r'Distance from Bunch Head S (m)', r'Scaled Wake Potential W (S)'],
                                     3: [r'Real Part of Longitudinal Impedance', r'Frequency f (GHz)', r'Real $Z_L$ ($k\Omega$) '],
                                     4: [r'Imaginary Part of Longitudinal Impedance', r'Frequency f (GHz)', r'Imag $Z_L$ ($k\Omega$) '],
                                     5: [r'Frequency Spectrum of Loss Factor', 'Frequency f (GHz)', r'$dk_L/df$ (V/pC/GHz) '],
                                     6: [r'Loss Factor Spectrum Integrated up to F', 'Frequency f (GHz)', r'$k_L$(F) (V/pC)'],
                                     7: [r'Real Part of Long. + Log Impedance', 'Frequency f (GHz)', r'Real $Z_{L+LOG}$ ($k\Omega$)'],
                                     8: [r'Imaginary Part of Long. + Log Impedance', 'Frequency f (GHz)', r'Imag $Z_{L+LOG}$ ($k\Omega$)'],
                                     9: [r'Spectrum of Long. + Log Loss Factor', 'Frequency f (GHz)', r'$dk_L/df$ ($V/pC/GHz$)'],
                                     10: [r'Long. + Log Factor Integrated up to F', 'Frequency f (GHz)', r'$k_L$(F) (V/pC)']
                                     }
        else:
            self.plot_decorations = {0: [r'Cavity Shape Input', r'Z-axis (m)', r'R-axis (m)'],
                                     1: [r'Cavity Shape Used', r'Z-axis (m)', r'R-axis (m)'],
                                     2: [r'Wake Potentials', r'Distance from Bunch Head S (m)', r'Scaled Wake Potential W (S)'],
                                     3: [r'Real Part of Azimuthal Impedance', r'Frequency f (GHz)', r'Real $Z_A$ ($k\Omega/m$) '],
                                     4: [r'Imaginary Part of Azimuthal Impedance', r'Frequency f (GHz)', r'Imag $Z_A$ ($k\Omega/m$) '],
                                     5: [r'Real Part of Transverse Impedance', r'Frequency f (GHz)', r'Real $Z_T$ ($k\Omega/m$) '],
                                     6: [r'Imaginary Part of Transverse Impedance', r'Frequency f (GHz)', r'Imag $Z_T$ ($k\Omega/m$) '],
                                     7: [r'Real Part of Longitudinal Impedance', r'Frequency f (GHz)', r'Real $Z_L$ ($k\Omega/m^2$) '],
                                     8: [r'Imaginary Part of Longitudinal Impedance', r'Frequency f (GHz)', r'Imag $Z_L$ ($k\Omega/m^2$) ']}

        return self.plot_decorations

    def _nparray_to_list(self, x):
        result = []
        for array in x:
            result.extend(array.tolist())

        return result

    def get_peaks(self, threshold):
        # Develop later
        peaks, _ = sps.find_peaks(self.y, height=threshold)
        xp, yp = np.array(self.x)[peaks], np.array(self.y)[peaks]

        return xp, yp
        # self.plot_beam_pipe_cutoffs()

    def _get_bands(self):
        try:
            # get band ranges
            group_dict = {}
            gr_n = 0
            xl, yl = [], []
            n = 0
            for i, val in enumerate(self.x_peaks):
                if i == 0:
                    xl.append(val)
                    yl.append(self.y_peaks[i])
                else:
                    # if val - xl[n] < 0.2:
                    if val - xl[n] < 0.1:
                        xl.append(val)
                        yl.append(self.y_peaks[i])
                        n += 1
                    else:
                        # save group
                        group_dict[gr_n] = {'xb': xl, 'yb': yl}

                        # next group
                        gr_n += 1

                        # reset lists
                        xl, yl = [], []
                        n = 0

                        # append current result
                        xl.append(val)
                        yl.append(self.y[i])

                    group_dict[gr_n] = {'xb': xl, 'yb': yl}  # last group

            # print_(group_dict)
            return group_dict
        except Exception as e:
            print_(f"Bands cannot be gotten without getting the peaks first. {e}")

    def _interp_fl(self):
        fl_list = []
        wx = [0, 0]
        wz = [1, 0]
        for f0, z in zip(self.x_peaks, self.y_peaks):
            z_0707 = 0.707 * z
            # print_(z_0707)

            l = 0
            for i, h in enumerate(self.y):
                if (l < z_0707 < h):

                    if (f0-0.01 < self.x[i] < f0):
                        wz = [l, h]
                        wx = [self.x[i - 1], self.x[i]]
                        break
                l = h

            # print_(f'\t{wx, wz}')

            # interpolate
            fl = wx[0] + (wx[1]-wx[0])*(z_0707-wz[0])/(wz[1]-wz[0])
            fl_list.append(fl)

        return fl_list

    def _interp_fr(self):
        fr_list = []
        wx = [0, 0]
        wz = [1, 0]
        for f0, z in zip(self.x_peaks, self.y_peaks):

            z_0707 = 0.707 * z
            # print_(z_0707)

            l = 0
            for i, h in enumerate(self.y):
                if l > z_0707 > h:
                    if (f0 < self.x[i] < f0 + 0.05):
                        wz = [l, h]
                        wx = [self.x[i - 1], self.x[i]]
                        break
                l = h

            # f0_old = f0
            # print_(f'\t{wx, wz}')

            # interpolate
            fr = wx[0] + (wx[1]-wx[0])*(z_0707-wz[0])/(wz[1]-wz[0])
            fr_list.append(fr)

        return fr_list

    def _get_Q(self):
        Q = [f0/bw for f0, bw in zip(self.x_peaks, self.BW)]
        return Q

    def _get_RQ(self):
        if self.MROT == 0:
            R_Q = [2*zz*1e3/q for zz, q in zip(self.y_peaks, self.Q)]
            # print_(f'R/Q: {R_Q}')
        else:
            c = 299792458 # m / s
            k = [2*np.pi*f0*1e9/c for f0 in self.x_peaks]
            R_Q_m = [2*zz*1e3/q for zz, q in zip(self.y_peaks, self.Q)]
            R_Q = [r_q_m/kk for r_q_m, kk in zip(R_Q_m, k)]
            # print_(f'R/Q: {R_Q}')

        return R_Q

    def _get_BW(self):
        BW = [a - b for a, b in zip(self.fr_list, self.fl_list)]
        return BW

    def plot_bands(self):
        fig, axs = plt.subplots(2, sharex=True)
        fig.suptitle('Vertically stacked subplots')
        axs[1].bar(self.x_peaks, self.R_Q, width=0.01, alpha=1)


        axs[0].plot(self.x, self.y)
        axs[0].scatter(self.x_peaks, self.y_peaks)
        color = ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c', '#ed718b', '#2ccf8b', '#fded21', '#8248d2', '#4c3a27', '#000000',
                 '#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c', '#ed718b', '#2ccf8b', '#fded21', '#8248d2', '#4c3a27', '#000000']

        for key, val in self.bands.items():
            tol = 0.05
            lb = val['xb'][0] - tol
            rb = val['xb'][-1] + tol
            axs[0].axvspan(lb, rb, 0, 1, color=color[key], alpha=0.3)
            axs[1].axvspan(lb, rb, 0, 1, color=color[key], alpha=0.3)

        plt.show()


class ABCIDataExtraction:
    def __init__(self):
        pass

    def multiple_folders_data(self, shape_space, abci_data_dir, request, save_excel, mon_interval=None, dip_interval=None, parallel=False):
        # process interval
        # Zmax
        if mon_interval is None:
            mon_interval = [[0.0, 2e10]]
        else:
            mon_interval = self.process_interval(mon_interval)

        if dip_interval is None:
            dip_interval = [[0.0, 2e10]]
        else:
            dip_interval = self.process_interval(dip_interval)

        reply = "Yes"
        def button_clicked(i):
            return i.text()

        if os.path.exists(f'{save_excel}.xlsx') and not parallel:
            msg = QMessageBox()
            msg.setWindowTitle("File Exist")
            msg.setText(f"Hey Chief, seems you've already processed the data for this folder. Do you want to overwrite it?")
            msg.setIcon(QMessageBox.Question)
            msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            msg.setDefaultButton(QMessageBox.No)

            msg.buttonClicked.connect(button_clicked)
            x = msg.exec_()

            if x == msg.Yes:
                reply = 'Yes'
            if x == msg.No:
                reply = 'No'

        if reply == "Yes":
            d = shape_space

            # data arrays initialization
            A, B, a, b, Ri, L, Req, alpha = [], [], [], [], [], [], [], []
            k_loss_array_transverse = []
            k_loss_array_longitudinal = []
            k_loss_M0 = []
            key_list = []

            # create list to hold Z
            Zmax_mon_list = []
            Zmax_dip_list = []
            for i in range(len(mon_interval)):
                Zmax_mon_list.append([])

            for i in range(len(dip_interval)):
                Zmax_dip_list.append([])

            def append_geom_parameters(values):
                A.append(values[0])
                B.append(values[1])
                a.append(values[2])
                b.append(values[3])
                Ri.append(values[4])
                L.append(values[5])
                Req.append(values[6])
                alpha.append(values[7])

            def calc_k_loss():
                for key, value in d.items():
                    print(f"Processing for Cavity {key}")
                    abci_data_long = ABCIData(abci_data_dir, key, 0)
                    abci_data_trans = ABCIData(abci_data_dir, key, 1)

                    # trans
                    x, y, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
                    k_loss_trans = abci_data_trans.loss_factor['Transverse']

                    if math.isnan(k_loss_trans):
                        print_(f"Encountered an exception: Check shape {key}")
                        continue

                    # long
                    x, y, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
                    abci_data_long.get_data('Loss Factor Spectrum Integrated up to F')

                    k_M0 = abci_data_long.y_peaks[0]
                    k_loss_long = abs(abci_data_long.loss_factor['Longitudinal'])
                    k_loss_HOM = k_loss_long - k_M0

                    # append only after successful run
                    append_geom_parameters(value)
                    key_list.append(key)
                    k_loss_M0.append(k_M0)
                    k_loss_array_longitudinal.append(k_loss_HOM)
                    k_loss_array_transverse.append(k_loss_trans)

            def get_Zmax(mon_interval=None, dip_interval=None):
                if mon_interval is None:
                    mon_interval = [0.0, 2e10]
                if dip_interval is None:
                    dip_interval = [0.0, 2e10]

                for key, value in d.keys():
                    print(f"Processing for Cavity {key}")
                    abci_data_mon = ABCIData(abci_data_dir, key, 0)
                    abci_data_dip = ABCIData(abci_data_dir, key, 1)

                    # get longitudinal and transverse impedance plot data
                    xr_mon, yr_mon, _ = abci_data_mon.get_data('Real Part of Longitudinal Impedance')
                    xi_mon, yi_mon, _ = abci_data_mon.get_data('Imaginary Part of Longitudinal Impedance')

                    xr_dip, yr_dip, _ = abci_data_dip.get_data('Real Part of Transverse Impedance')
                    xi_dip, yi_dip, _ = abci_data_dip.get_data('Imaginary Part of Transverse Impedance')

                    # Zmax
                    if mon_interval is None:
                        mon_interval = [[0.0, 10]]
                    if dip_interval is None:
                        dip_interval = [[0.0, 10]]

                    # calculate magnitude
                    ymag_mon = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_mon, yi_mon)]
                    ymag_dip = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_dip, yi_dip)]

                    # get peaks
                    peaks_mon, _ = sps.find_peaks(ymag_mon, height=0)
                    xp_mon, yp_mon = np.array(xr_mon)[peaks_mon], np.array(ymag_mon)[peaks_mon]

                    peaks_dip, _ = sps.find_peaks(ymag_dip, height=0)
                    xp_dip, yp_dip = np.array(xr_dip)[peaks_dip], np.array(ymag_dip)[peaks_dip]

                    for i, z_bound in enumerate(mon_interval):
                        # get mask
                        msk_mon = [(z_bound[0] < x < z_bound[1]) for x in xp_mon]

                        if len(yp_mon[msk_mon]) != 0:
                            Zmax_mon = max(yp_mon[msk_mon])

                            Zmax_mon_list[i].append(Zmax_mon)
                        elif len(yp_mon) != 0:
                            Zmax_mon_list[i].append(0)
                        else:
                            continue

                    for i, z_bound in enumerate(dip_interval):
                        # get mask
                        msk_dip = [(z_bound[0] < x < z_bound[1]) for x in xp_dip]

                        if len(yp_dip[msk_dip]) != 0:
                            Zmax_dip = max(yp_dip[msk_dip])

                            Zmax_dip_list[i].append(Zmax_dip)
                        elif len(yp_dip) != 0:
                            Zmax_dip_list[i].append(0)
                        else:
                            continue

                    append_geom_parameters(value["IC"])

            def all(mon_interval, dip_interval):
                for key, value in d.items():
                    print(f"Processing for Cavity {key}")
                    abci_data_long = ABCIData(abci_data_dir, key, 0)
                    abci_data_trans = ABCIData(abci_data_dir, key, 1)

                    # get longitudinal and transverse impedance plot data
                    xr_mon, yr_mon, _ = abci_data_long.get_data('Real Part of Longitudinal Impedance')
                    xi_mon, yi_mon, _ = abci_data_long.get_data('Imaginary Part of Longitudinal Impedance')

                    xr_dip, yr_dip, _ = abci_data_trans.get_data('Real Part of Transverse Impedance')
                    xi_dip, yi_dip, _ = abci_data_trans.get_data('Imaginary Part of Transverse Impedance')

                    # loss factors
                    # trans
                    k_loss_trans = abci_data_trans.loss_factor['Transverse']

                    if math.isnan(k_loss_trans):
                        print_(f"Encountered an exception: Check shape {key}")
                        continue

                    # long
                    abci_data_long.get_data('Loss Factor Spectrum Integrated upto F')

                    k_M0 = abci_data_long.y_peaks[0]
                    k_loss_long = abs(abci_data_long.loss_factor['Longitudinal'])
                    k_loss_HOM = k_loss_long - k_M0

                    # calculate magnitude
                    ymag_mon = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_mon, yi_mon)]
                    ymag_dip = [(a ** 2 + b ** 2) ** 0.5 for a, b in zip(yr_dip, yi_dip)]

                    # get peaks
                    peaks_mon, _ = sps.find_peaks(ymag_mon, height=0)
                    xp_mon, yp_mon = np.array(xr_mon)[peaks_mon], np.array(ymag_mon)[peaks_mon]

                    peaks_dip, _ = sps.find_peaks(ymag_dip, height=0)
                    xp_dip, yp_dip = np.array(xr_dip)[peaks_dip], np.array(ymag_dip)[peaks_dip]

                    for i, z_bound in enumerate(mon_interval):
                        # get mask
                        msk_mon = [(z_bound[0] < x < z_bound[1]) for x in xp_mon]

                        if len(yp_mon[msk_mon]) != 0:
                            Zmax_mon = max(yp_mon[msk_mon])

                            Zmax_mon_list[i].append(Zmax_mon)
                        elif len(yp_mon) != 0:
                            Zmax_mon_list[i].append(0)
                        else:
                            continue

                    for i, z_bound in enumerate(dip_interval):
                        # get mask
                        msk_dip = [(z_bound[0] < x < z_bound[1]) for x in xp_dip]

                        if len(yp_dip[msk_dip]) != 0:
                            Zmax_dip = max(yp_dip[msk_dip])

                            Zmax_dip_list[i].append(Zmax_dip)
                        elif len(yp_dip) != 0:
                            Zmax_dip_list[i].append(0)
                        else:
                            continue

                    # append only after successful run
                    append_geom_parameters(value["IC"])
                    key_list.append(key)

                    k_loss_M0.append(k_M0)
                    k_loss_array_longitudinal.append(k_loss_HOM)
                    k_loss_array_transverse.append(k_loss_trans)

            if request.lower() == 'k_loss':
                calc_k_loss()
                # save excel
                try:
                    data = {'key': key_list, 'A': A, 'B': B, 'a': a, 'b': b, 'Ri': Ri, 'L': L, 'Req': Req, "alpha": alpha,
                            'k_loss_M0': k_loss_M0, 'k_loss_long': k_loss_array_longitudinal, 'k_loss_trans': k_loss_array_transverse}

                    df = pd.DataFrame.from_dict(data)
                    df.to_excel(f'{save_excel}.xlsx', index=False)
                except Exception as e:
                    print("Oops! Encountered some error trying to save file: ", e)

            if request.lower() == 'zmax':
                get_Zmax(mon_interval, dip_interval)
                print(len(Zmax_mon_list), len(Zmax_dip_list), len(A))

                # save excel
                try:
                    data = {'key': key_list, 'A': A, 'B': B, 'a': a, 'b': b, 'Ri': Ri, 'L': L, 'Req': Req, "alpha": alpha}

                    # add impedance lists
                    for i, msk in enumerate(mon_interval):
                        data[f'Z_long[max({msk[0]}<f<{msk[1]})]'] = Zmax_mon_list[i]

                    for i, msk in enumerate(dip_interval):
                        data[f'Z_trans[max({msk[0]}<f<{msk[1]})]'] = Zmax_dip_list[i]

                    df = pd.DataFrame.from_dict(data)
                    df.to_excel(f'{save_excel}.xlsx', index=False)
                except Exception as e:
                    print("Oops! Encountered some error trying to save file: ", e)

            if request.lower() == 'all':
                all(mon_interval, dip_interval)

                # save excel
                try:
                    data = {'key': key_list, 'A': A, 'B': B, 'a': a, 'b': b, 'Ri': Ri, 'L': L, 'Req': Req, "alpha": alpha,
                            'k_loss_M0': k_loss_M0, 'k_loss_long': k_loss_array_longitudinal, 'k_loss_trans': k_loss_array_transverse}

                    # add impedance lists
                    for i, msk in enumerate(mon_interval):
                        data[f'Z_long[max({msk[0]}<f<{msk[1]})]'] = Zmax_mon_list[i]

                    for i, msk in enumerate(dip_interval):
                        data[f'Z_trans[max({msk[0]}<f<{msk[1]})]'] = Zmax_dip_list[i]

                    df = pd.DataFrame.from_dict(data)
                    df.to_excel(f'{save_excel}.xlsx', index=False)

                except Exception as e:
                    print("Oops! Encountered some error trying to save file: ", e)

    def multiple_folders_data_parallel(self, shape_space, abci_data_folder, proc_count, request, save_excel, temp_folder, mon_interval=None, dip_interval=None):

        # create temporary folder
        if os.path.exists(fr"{temp_folder}"):
            pass
        else:
            os.mkdir(fr"{temp_folder}")

        processes = []

        keys = list(shape_space.keys())
        shape_space_len = len(keys)
        share = round(shape_space_len / proc_count)

        for p in range(proc_count):
            # create temporary shape spaces
            if p < proc_count - 1:
                proc_keys_list = keys[p * share:p * share + share]
            else:
                proc_keys_list = keys[p * share:]

            print(proc_keys_list)
            processor_shape_space = {}

            for key, val in shape_space.items():
                if key in proc_keys_list:
                    processor_shape_space[key] = val

            service = mp.Process(target=self.multiple_folders_data,
                                 args=(processor_shape_space, abci_data_folder, request),
                                 # keys[-1]]
                                 kwargs={'save_excel': f'{temp_folder}\Proc_{p}', 'mon_interval': mon_interval,
                                         'dip_interval': dip_interval, 'parallel': True})
            service.start()
            processes.append(service)

        for p in processes:
            p.join()

        # join temporary files and delete
        self.join_excel('Proc', proc_count, save_excel, temp_folder)

        # delete temporary folder
        shutil.rmtree(temp_folder)

    def process_interval(self, interval_list):
        interval = []
        for i in range(len(interval_list)-1):
            interval.append([interval_list[i], interval_list[i+1]])

        return interval

    def join_excel(self, generic_name, proc_count, save_excel, temp_folder):
        df = fr.excel_reader(fr'{temp_folder}\{generic_name}_{0}.xlsx')['Sheet1']

        for p in range(1, proc_count):
            d = fr.excel_reader(fr'{temp_folder}\{generic_name}_{p}.xlsx')
            d = d['Sheet1']
            df = pd.merge(df, d, how='outer')

        try:
            df.to_excel(f'{save_excel}.xlsx', index=False)
        except Exception as e:
            print("Oops! Encountered some error trying to save file: ", e)


if __name__ == '__main__':
    directory = r"D:\Dropbox\Projects\NewFolder\SimulationData\ABCI"
    fid = "Cavity0"  # folder name
    MROT = "1"  # 0 for monopole, 1 for dipole

    # create ABCIData object
    abci_data = ABCIData(directory, fid, MROT)

    # call get_data method and provide key
    keys = {0: r'Cavity Shape Input',
            1: r'Cavity Shape Used',
            2: r'Wake Potentials',
            3: r'Real Part of Longitudinal Impedance',
            4: r'Imaginary Part of Longitudinal Impedance',
            5: 'Frequency Spectrum of Loss Factor',
            6: r'Loss Factor Spectrum Integrated upto F',
            7: r'Real Part of Long. + Log Impedance',
            8: r'Imaginary Part of Long. + Log Impedance',
            9: r'Spectrum of Long. + Log Loss Factor',
            10: r'Long. + Log Factor Integrated up to F',
            11: r'Real Part of Azimuthal Impedance',
            12: r'Imaginary Part of Azimuthal Impedance',
            13: r'Real Part of Transverse Impedance',
            14: r'Imaginary Part of Transverse Impedance'}

    x, y, _ = abci_data.get_data(key=3)  # For the key, either the title of the plot can be given as input or the index
    print(x)
    print(y)

    import matplotlib.pyplot as plt
    plt.plot(x, y)
    plt.show()

