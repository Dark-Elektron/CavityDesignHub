import json
import os
import shutil

import scipy
from PyQt5.QtWidgets import QMessageBox
import multiprocessing as mp

from icecream import ic

from utils.file_reader import FileReader
import numpy as np
import pandas as pd

from utils.shared_functions import perform_geometry_checks

fr = FileReader()


class SLANSData:
    ##### COMPLETE CODE
    def __init__(self, dirc, fid, bc):
        self.mode_count = None
        self.title_dict = {}
        self.data_dict = {}
        self.dir = dirc
        self.fid = fid
        self.bc = bc

        self.path = dirc
        # print(self.path)

    def get_0D_plot_data(self):
        # print("\t\tSTARTED WITH GETTING 0D PLOT DATA")
        filename = fr'{self.path}\{self.fid}\cavity_{self.bc}.svl'
        # print(filename)

        data_dict = fr.svl_reader(filename)
        self.mode_count = len(data_dict['FREQUENCY'])
        return data_dict

    def get_1D_plot_data(self):
        # print("\t\tSTARTED WITH GETTING 1D PLOT DATA")
        data_dict = {}
        # SET MODE COUNT
        self.get_0D_plot_data()

        # print(self.mode_count)
        for mode in range(1, self.mode_count + 1):
            data_dict[mode] = {}
            data_dict[mode] = {}

            data_dict[mode]["Axis Fields"] = self.get_axis_field_data(mode)
            data_dict[mode]["Surface Fields"] = self.get_surface_field_data(mode)

        # print("\t\tDONE WITH GETTING 1D PLOT DATA")
        return data_dict

    def get_1D_plot_title(self):
        pass

    def get_2D_plot_data(self):
        print("\t\tSTARTED WITH GETTING 2D PLOT DATA")
        data_dict = {}
        # SET MODE COUNT
        self.get_0D_plot_data()

        print(self.mode_count)
        # FOR NOW, ONLY ONE MODE RESULT IS EXTRACTABLE, CHANGE LATER
        for mode in range(1, 2):
            data_dict[mode] = {}
            data_dict[mode] = {}

            fields_on_plane_data = {}
            Ez, Er, Emag, Hi = [], [], [], []
            print(fr"{self.path}/Cavity{self.fid}/cavity_{self.bc}.pam")
            with open(fr"{self.path}/Cavity{self.fid}/cavity_{self.bc}.pam", 'r') as f:
                # print(f.read())
                for l in f:
                    # print([float(x) for x in l.strip().split(' ') if x !=''])
                    ll = [float(x) for x in l.strip().split(' ') if x != '']

                    Ez.append(ll[0])
                    Er.append(ll[1])
                    Emag.append(ll[2])
                    #### line 2 has just three entries
                    if len(ll) < 4:
                        Hi.append(0)
                    else:
                        Hi.append(ll[3])

            #### GET R AND Z LIMITS AND INTERVALS
            print("\t\tGetting r and z")
            Z = Er[0]
            R = Er[1]
            nZ = Emag[0]
            nR = Emag[1]

            Ez = np.array(Ez)[2:]
            Er = np.array(Er)[2:]
            Emag = np.array(Emag)[2:]
            Hi = np.array(Hi)[2:]
            # print(Ez)

            #### CREATE MESH GRID
            print("Crete mesh grid")
            R = np.linspace(0, R, int(nR) + 1)
            Z = np.linspace(0, Z, int(nZ) + 1)
            Rmsh, Zmsh = np.meshgrid(Z, R)

            #### DIRECTION GRID AND DATA
            print("Direction Grid node_editor")
            Er_mat = Er.reshape(len(R), len(Z))
            Ez_mat = Ez.reshape(len(R), len(Z))
            Emag_mat = Emag.reshape(len(R), len(Z))
            Hi_mat = Hi.reshape(len(R), len(Z))

            ### MASKING
            print("Masking")
            Emag_mat_masked = np.ma.array(Emag_mat)

            # mask values below a certain threshold
            Ez_mat_masked = np.ma.masked_where(Emag_mat_masked == 0, Ez_mat)
            Er_mat_masked = np.ma.masked_where(Emag_mat_masked == 0, Er_mat)
            Emag_mat_masked = np.ma.masked_where(Emag_mat_masked == 0, Emag_mat)
            Hi_mat_masked = np.ma.masked_where(Emag_mat_masked == 0, Hi_mat)

            Ez_dict = {'R': Rmsh, 'Z': Zmsh, "Ez": Ez_mat_masked}
            Er_dict = {'R': Rmsh, 'Z': Zmsh, "Er": Er_mat_masked}
            Emag_dict = {'R': Rmsh, 'Z': Zmsh, "Emag": Emag_mat_masked}
            Hi_dict = {'R': Rmsh, 'Z': Zmsh, "Hi": Hi_mat_masked}

            fields_on_plane_data['Ez'] = Ez_dict
            fields_on_plane_data['Er'] = Er_dict
            fields_on_plane_data['Emag'] = Emag_dict
            fields_on_plane_data['Hi'] = Hi_dict

            data_dict[mode] = fields_on_plane_data

        print("\t\tDONE WITH GETTING 2D PLOT DATA")
        return data_dict
        # plt.set_cmap("jet")
        # plt.contourf(Rmsh, Zmsh, Emag_mat_masked)
        # plt.quiver(Rmsh, Zmsh, Ez_mat, Er_mat, Emag_mat_masked, scale_units='width', scale=5)
        # plt.show()

    def get_2D_plot_title(self):
        pass

    def process_line(self, line):
        line = line.strip().split(' ')
        res = 0
        for val in line:
            try:
                res = float(val)
                break
            except:
                continue
        return res

    def get_surface_field_data(self, mode):
        # SURFACE FIELDS
        y = []
        surface_field_dict = {}
        with open(fr"{self.path}\Cavity{self.fid}\cavity_{self.bc}_{mode}.sf", 'r') as f:
            for ll in f.readlines():
                ll = ll.strip()
                y.append(abs(float(ll)))

        surface_field_dict['Es'] = y

        return surface_field_dict

    def get_axis_field_data(self, mode):
        axis_field_data = {}

        x, y = [], []
        path = os.path.join(os.getcwd(), fr"{self.path}\Cavity{self.fid}\cavity_{self.bc}_{mode}.af")
        with open(path, 'r') as f:
            for ll in f.readlines():
                ll = ll.strip()
                x.append(float(ll.split(' ')[0]))
                y.append(float(ll.split(' ')[1]))

        # get avg x
        # avg_x = sum(x) / len(x)
        # x_shift = [t - avg_x for t in x]

        y_abs = [abs(e) for e in y]

        # RETURN ABSOLUTE FIELD VALUE
        axis_field_data['x'], axis_field_data['y'], axis_field_data['y_abs'] = x, y, y_abs
        return axis_field_data
        # y_norm = [abs(t) / max(y_abs) for t in y_abs]

        # ax.plot(x_shift, y_abs)

    def normal_interp(self, x, y, a, xi, yi):
        rbf = scipy.interpolate.Rbf(x, y, a)
        ai = rbf(xi, yi)
        return ai

    def get_mode_count(self):
        self.get_0D_plot_data()
        return self.mode_count


class SLANSDataExtraction:
    def __init__(self):
        pass

    def multiple_folders_data(self, dirs, slans_data_folder, mode, bc, request, save_excel, parallel=False):
        reply = "Yes"

        def button_clicked(i):
            return i.text()

        if os.path.exists(f'{save_excel}.xlsx') and not parallel:
            msg = QMessageBox()
            msg.setWindowTitle("File Exist")
            msg.setText(
                f"Hey Chief, seems you've already processed the data for this folder. Do you want to overwrite it?")
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
            # build shape space
            df = pd.DataFrame(columns=['key',
                                       'A_i', 'B_i', 'a_i', 'b_i', 'Ri_i', 'L_i', 'Req', 'alpha_i',
                                       "Req [mm]", "Normalization Length [mm]", "N cells", "freq [MHz]", "Q []",
                                       "E [MV/m]", "Vacc [MV]", "Eacc [MV/m]", "Epk [MV/m]", "Hpk [A/m]", "Bpk [mT]",
                                       "kcc [%]", "ff [%]", "Rsh [MOhm]", "R/Q [Ohm]",
                                       "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "G [Ohm]", "GR/Q [Ohm^2]"])

            for key in dirs:
                if os.path.exists(f"{slans_data_folder}/{key}/qois.json") \
                        and os.path.exists(f"{slans_data_folder}/{key}/geometric_parameters.json"):
                    geometric_parameter = {}
                    qois = {}

                    with open(f"{slans_data_folder}/{key}/geometric_parameters.json") as json_file:
                        geometric_parameter.update(json.load(json_file))

                    with open(f"{slans_data_folder}/{key}/qois.json") as json_file:
                        qois.update(json.load(json_file))

                    # check frequency
                    if not 801 <= qois['freq [MHz]'] <= 802:
                        continue

                    # check alpha
                    if not perform_geometry_checks(geometric_parameter['IC'], geometric_parameter['OC']):
                        continue


                    df.loc[key, ['key']] = key
                    df.loc[key, ['A_i', 'B_i', 'a_i', 'b_i', 'Ri_i', 'L_i', 'Req', 'alpha_i']] = geometric_parameter[
                        'IC']

                    # check if end cell is not the same as mid cell
                    if geometric_parameter['IC'] != geometric_parameter['OC']:
                        df.loc[key, ['Aeq_ol', 'Beq_ol', 'ai_ol', 'bi_ol', 'Ri_ol', 'L_ol', 'Req_ol', 'alpha_ol']] = \
                        geometric_parameter['OC']

                        if 'OC_R' in list(geometric_parameter.keys()):
                            if geometric_parameter['OC'] != geometric_parameter['OC_R']:
                                df.iloc[
                                    key, ['Aeq_or', 'Beq_or', 'ai_or', 'bi_or', 'Ri_or', 'L_or', 'Req_or',
                                          'alpha_or']] = \
                                    geometric_parameter['OC']
                    df.loc[key, ["Req [mm]", "Normalization Length [mm]", "N cells", "freq [MHz]", "Q []",
                                 "E [MV/m]", "Vacc [MV]", "Eacc [MV/m]", "Epk [MV/m]", "Hpk [A/m]", "Bpk [mT]",
                                 "kcc [%]", "ff [%]", "Rsh [MOhm]", "R/Q [Ohm]",
                                 "Epk/Eacc []", "Bpk/Eacc [mT/MV/m]", "G [Ohm]", "GR/Q [Ohm^2]"]] = list(qois.values())

                # print(key)
            ic(df)
            df = df.sort_values(by='key')
            # df = pd.DataFrame.from_dict(data)
            df.to_excel(f'{save_excel}.xlsx', index=False)


            # d = shape_space
            #
            # # data arrays initialization
            # A, B, a, b, Ri, L, Req, alpha = [], [], [], [], [], [], [], []
            # key_list = []
            #
            # # data arrays initialization
            # Freq_arr = []
            # E_stored_arr = []
            # Rsh_arr = []
            # Q_arr = []
            # Epk_arr = []
            # Hpk_arr = []
            # Eacc_arr = []
            # Rsh_Q_arr = []
            # Epk_Eacc_arr = []
            # Bpk_Eacc_arr = []
            # kcc_arr = []
            #
            # def append_geom_parameters(values):
            #     A.append(values[0])
            #     B.append(values[1])
            #     a.append(values[2])
            #     b.append(values[3])
            #     Ri.append(values[4])
            #     L.append(values[5])
            #     Req.append(values[6])
            #     alpha.append(values[7])
            #
            # def get_0D_plot_data():
            #
            #     # directory_list = [os.path.join(slans_data_dir, name) for name in os.listdir(slans_data_dir) if
            #     #                   os.path.isdir(os.path.join(slans_data_dir, name))]
            #     # print(directory_list)
            #     for key, value in d.items():
            #         # if os.path.exists(f"{dir_}\geometric_parameters.json"):
            #         #     with open(f"{dir_}\geometric_parameters.json", 'r') as f:
            #         #         d = json.load(f)
            #         # else:
            #         #     continue
            #         # key = dir_[dir_.index("Cavity"):]
            #         slans_data = SLANSData(slans_data_dir, f"{key}", bc)
            #
            #         try:
            #             d_0d = slans_data.get_0D_plot_data()
            #
            #             Freq = d_0d['FREQUENCY'][mode - 1]
            #             E_stored = d_0d['STORED ENERGY'][mode - 1]
            #             Rsh = d_0d['SHUNT IMPEDANCE'][mode - 1]  # MOhm
            #             Q = d_0d['QUALITY FACTOR'][mode - 1]
            #             Epk = d_0d['MAXIMUM ELEC. FIELD'][mode - 1]  # MV/m
            #             Hpk = d_0d['MAXIMUM MAG. FIELD'][mode - 1]  # A/m
            #             # Vacc = d_0d['ACCELERATION'][mode-1]
            #             Eavg = d_0d['AVERAGE E.FIELD ON AXIS'][mode - 1]  # MV/m
            #             Rsh_Q = d_0d['EFFECTIVE IMPEDANCE'][mode - 1]  # Ohm
            #             # Epk_Eacc = d_0d['KM (Emax/Accel.rate)'][mode-1]  #
            #             # Bpk_Eacc = d_0d['KH (Hmax*Z0/Accel.rate)'][mode-1]  # mT/(MV/m)
            #             try:
            #                 L = 4*value["IC"][5]  # normalize to length of cell
            #             except KeyError:
            #                 L = 4*value['L']
            #
            #             Vacc = np.sqrt(2*Rsh_Q * E_stored * 2*np.pi * Freq * 1e6) * 1e-6  # factor of 2, remember circuit and accelerator definition
            #             # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
            #             # Eacc = Vacc / (374 * 1e-3)  # for 1 cell factor of 2, remember circuit and accelerator definition
            #             Eacc = Vacc / (L * 1e-3)  # for 1 cell factor of 2, remember circuit and accelerator definition
            #             Epk_Eacc = Epk / Eacc
            #             Bpk_Eacc = (Hpk * 4 * np.pi * 1e-7) * 1e3 / Eacc
            #
            #             Freq_arr.append(Freq)
            #             E_stored_arr.append(E_stored)
            #             Rsh_arr.append(2 * Rsh)  # corrected to accelerator definition
            #             Q_arr.append(Q)
            #             Epk_arr.append(Epk)
            #             Hpk_arr.append(Hpk)
            #             Eacc_arr.append(Eacc)
            #             Epk_Eacc_arr.append(Epk_Eacc)
            #             Bpk_Eacc_arr.append(Bpk_Eacc)
            #             Rsh_Q_arr.append(2 * Rsh_Q)  # corrected to accelerator definition
            #
            #             f_0, f_pi = d_0d['FREQUENCY'][0], d_0d['FREQUENCY'][mode - 1]
            #             try:
            #                 kcc = 2 * (f_pi - f_0) / (f_pi + f_0)
            #             except ZeroDivisionError:
            #                 kcc = 0
            #
            #             kcc_arr.append(kcc * 100)
            #             try:
            #                 append_geom_parameters(value["IC"])
            #             except KeyError:
            #                 append_geom_parameters(list(value.values()))
            #
            #             key_list.append(key)
            #
            #         except FileNotFoundError:
            #             pass
            #
            # # sheets = ['Sheet1', 'Sheet2']
            # def get_1D_plot_data():
            #     new_save = True
            #     for key, value in d.items():
            #
            #         slans_data = SLANSData(slans_data_dir, key, bc)
            #         d_1d = slans_data.get_1D_plot_data()
            #
            #         dd = d_1d[mode]
            #         fields = list(dd.keys())
            #         # print(sheets)
            #
            #         # save excel
            #         if save_excel:
            #             data = {'x': dd[fields[0]]['x'], 'y_abs': dd[fields[0]]['y_abs']}
            #             df_axis_field = pd.DataFrame.from_dict(data)
            #
            #             data = {'Es': dd[fields[1]]['Es']}
            #             df_surface_field = pd.DataFrame.from_dict(data)
            #
            #             if new_save:
            #                 with pd.ExcelWriter(f'{save_excel}_Axis_Fields.xlsx') as writer:
            #                     df_axis_field.to_excel(writer, sheet_name=f'{key}')
            #
            #                 with pd.ExcelWriter(f'{save_excel}_Surface_Fields.xlsx') as writer:
            #                     df_surface_field.to_excel(writer, sheet_name=f'{key}')
            #
            #                 new_save = False
            #             else:
            #                 with pd.ExcelWriter(f'{save_excel}_Axis_Fields.xlsx', engine="openpyxl",
            #                                     mode='a') as writer:
            #                     df_axis_field.to_excel(writer, sheet_name=f'{key}')
            #
            #                 with pd.ExcelWriter(f'{save_excel}_Surface_Fields.xlsx', engine="openpyxl",
            #                                     mode='a') as writer:
            #                     df_surface_field.to_excel(writer, sheet_name=f'{key}')
            #
            # print(request)
            #
            # if request == '0D Data':
            #     get_0D_plot_data()
            #     print(len(Epk_Eacc_arr), len(Bpk_Eacc_arr), len(Epk_arr), len(A), len(kcc_arr))
            #
            #     # save excel
            #     if save_excel:
            #         data = {'key': key_list, 'A': A, 'B': B, 'a': a, 'b': b, 'Ri': Ri, 'L': L, 'Req': Req,
            #                 "alpha": alpha,
            #                 f'Freq_{mode}': Freq_arr,
            #                 f'E_stored_{mode}': E_stored_arr, f'Rsh_{mode}': Rsh_arr, f'Q_{mode}': Q_arr,
            #                 f'Epk_{mode}': Epk_arr, f'Hpk_{mode}': Hpk_arr, f'Eacc_{mode}': Eacc_arr,
            #                 f'Rsh/Q_{mode}': Rsh_Q_arr, f'Epk/Eacc_{mode}': Epk_Eacc_arr,
            #                 f'Bpk/Eacc_{mode}': Bpk_Eacc_arr, f'kcc_{mode}': kcc_arr}
            #
            #         df = pd.DataFrame.from_dict(data)
            #         df.to_excel(f'{save_excel}.xlsx', index=False)
            #
            #         # # save to json
            #         # with open(f'{save_excel}.json', 'w') as f:
            #         #     json.dump(data, f, indent=4, separators=(',', ': '))
            #
            # if request == '1D Data':
            #     if not os.path.exists(f'{save_excel}_Axis_Fields.xlsx') and not os.path.exists(
            #             f'{save_excel}_Surface_Fields.xlsx'):
            #         get_1D_plot_data()
            #     else:
            #         print(f"Hey Chief, seems you've already processed the {save_excel} data for this folder. "
            #               "Delete the file to reprocess or rename save_excel argument")

    def multiple_folders_data_old(self, shape_space, slans_data_dir, mode, bc, request, save_excel, parallel=False):
        reply = "Yes"

        def button_clicked(i):
            return i.text()

        if os.path.exists(f'{save_excel}.xlsx') and not parallel:
            msg = QMessageBox()
            msg.setWindowTitle("File Exist")
            msg.setText(
                f"Hey Chief, seems you've already processed the data for this folder. Do you want to overwrite it?")
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
            key_list = []

            # data arrays initialization
            E_stored_arr = []
            Rsh_arr = []
            Q_arr = []
            Epk_arr = []
            Hpk_arr = []
            Eacc_arr = []
            Rsh_Q_arr = []
            Epk_Eacc_arr = []
            Bpk_Eacc_arr = []
            kcc_arr = []

            def append_geom_parameters(values):
                A.append(values[0])
                B.append(values[1])
                a.append(values[2])
                b.append(values[3])
                Ri.append(values[4])
                L.append(values[5])
                Req.append(values[6])
                alpha.append(values[7])

            def get_0D_plot_data():
                for key, value in d.items():
                    slans_data = SLANSData(slans_data_dir, key, bc)

                    d_0d = slans_data.get_0D_plot_data()

                    E_stored = d_0d['STORED ENERGY'][mode - 1]
                    Rsh = d_0d['SHUNT IMPEDANCE'][mode - 1]  # MOhm
                    Q = d_0d['QUALITY FACTOR'][mode - 1]
                    Epk = d_0d['MAXIMUM ELEC. FIELD'][mode - 1]  # MV/m
                    Hpk = d_0d['MAXIMUM MAG. FIELD'][mode - 1]  # A/m
                    # Vacc = d_0d['ACCELERATION'][mode-1]
                    Eavg = d_0d['AVERAGE E.FIELD ON AXIS'][mode - 1]  # MV/m
                    Rsh_Q = d_0d['EFFECTIVE IMPEDANCE'][mode - 1]  # Ohm
                    # Epk_Eacc = d_0d['KM (Emax/Accel.rate)'][mode-1]  #
                    # Bpk_Eacc = d_0d['KH (Hmax*Z0/Accel.rate)'][mode-1]  # mT/(MV/m)

                    Vacc = np.sqrt(
                        Rsh_Q * E_stored * np.pi * 400.79 * 1e6) * 1e-6  # factor of 2, remember circuit and accelerator definition
                    Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
                    Epk_Eacc = Epk / (Eacc)
                    Bpk_Eacc = (Hpk * 4 * np.pi * 1e-7) * 1e3 / (Eacc)

                    E_stored_arr.append(E_stored)
                    Rsh_arr.append(2 * Rsh)  # corrected to accelerator definition
                    Q_arr.append(Q)
                    Epk_arr.append(Epk)
                    Hpk_arr.append(Hpk)
                    Eacc_arr.append(Eacc)
                    Epk_Eacc_arr.append(Epk_Eacc)
                    Bpk_Eacc_arr.append(Bpk_Eacc)
                    Rsh_Q_arr.append(2 * Rsh_Q)  # corrected to accelerator definition

                    f_0, f_pi = d_0d['FREQUENCY'][0], d_0d['FREQUENCY'][mode - 1]
                    kcc = 2 * (f_pi - f_0) / (f_pi + f_0)
                    kcc_arr.append(kcc * 100)
                    append_geom_parameters(value["IC"])
                    key_list.append(key)

            # sheets = ['Sheet1', 'Sheet2']
            def get_1D_plot_data():
                new_save = True
                for key, value in d.items():

                    slans_data = SLANSData(slans_data_dir, key, bc)
                    d_1d = slans_data.get_1D_plot_data()

                    dd = d_1d[mode]
                    fields = list(dd.keys())
                    # print(sheets)

                    # save excel
                    if save_excel:
                        data = {'x': dd[fields[0]]['x'], 'y_abs': dd[fields[0]]['y_abs']}
                        df_axis_field = pd.DataFrame.from_dict(data)

                        data = {'Es': dd[fields[1]]['Es']}
                        df_surface_field = pd.DataFrame.from_dict(data)

                        if new_save:
                            with pd.ExcelWriter(f'{save_excel}_Axis_Fields.xlsx') as writer:
                                df_axis_field.to_excel(writer, sheet_name=f'{key}')

                            with pd.ExcelWriter(f'{save_excel}_Surface_Fields.xlsx') as writer:
                                df_surface_field.to_excel(writer, sheet_name=f'{key}')

                            new_save = False
                        else:
                            with pd.ExcelWriter(f'{save_excel}_Axis_Fields.xlsx', engine="openpyxl",
                                                mode='a') as writer:
                                df_axis_field.to_excel(writer, sheet_name=f'{key}')

                            with pd.ExcelWriter(f'{save_excel}_Surface_Fields.xlsx', engine="openpyxl",
                                                mode='a') as writer:
                                df_surface_field.to_excel(writer, sheet_name=f'{key}')

            print(request)

            if request == '0D node_editor':
                get_0D_plot_data()
                print(len(Epk_Eacc_arr), len(Bpk_Eacc_arr), len(Epk_arr), len(A), len(kcc_arr))

                # save excel
                if save_excel:
                    data = {'key': key_list, 'A': A, 'B': B, 'a': a, 'b': b, 'Ri': Ri, 'L': L, 'Req': Req,
                            "alpha": alpha,
                            f'E_stored_{mode}': E_stored_arr, f'Rsh_{mode}': Rsh_arr, f'Q_{mode}': Q_arr,
                            f'Epk_{mode}': Epk_arr, f'Hpk_{mode}': Hpk_arr, f'Eacc_{mode}': Eacc_arr,
                            f'Rsh/Q_{mode}': Rsh_Q_arr, f'Epk/Eacc_{mode}': Epk_Eacc_arr,
                            f'Bpk/Eacc_{mode}': Bpk_Eacc_arr, f'kcc_{mode}': kcc_arr}

                    df = pd.DataFrame.from_dict(data)
                    df.to_excel(f'{save_excel}.xlsx', index=False)

            if request == '1D node_editor':
                if not os.path.exists(f'{save_excel}_Axis_Fields.xlsx') and not os.path.exists(
                        f'{save_excel}_Surface_Fields.xlsx'):
                    get_1D_plot_data()
                else:
                    print(f"Hey Chief, seems you've already processed the {save_excel} data for this folder. "
                          "Delete the file to reprocess or rename save_excel argument")

    def multiple_folders_data_parallel(self, slans_data_folder, proc_count, mode, bc, request, save_excel,
                                       temp_folder):

        # bc = bc.replace('m', '3')
        # bc = bc.replace('e', '2')

        # create temporary folder
        if os.path.exists(fr"{temp_folder}"):
            pass
        else:
            os.mkdir(fr"{temp_folder}")

        processes = []

        # get list of directories in folder
        dirs = [x for x in os.listdir(slans_data_folder) if '_process_' not in x]
        # print(dirs)


        keys = dirs
        shape_space_len = len(keys)
        share = round(shape_space_len / proc_count)

        for p in range(proc_count):
            # create temporary shape spaces
            if p < proc_count - 1:
                proc_keys_list = keys[p * share:p * share + share]
            else:
                proc_keys_list = keys[p * share:]

            # print(proc_keys_list)
            # processor_shape_space = {}
            #
            # for key, val in shape_space.items():
            #     if key in proc_keys_list:
            #         processor_shape_space[key] = val

            service = mp.Process(target=self.multiple_folders_data,
                                 args=(proc_keys_list, slans_data_folder, mode, bc, request),
                                 kwargs={'save_excel': f'{temp_folder}\Proc_{p}', 'parallel': True})
            service.start()
            processes.append(service)

        for p in processes:
            p.join()

        # join temporary files and delete
        self.join_excel('Proc', proc_count, save_excel, temp_folder)
        print('Done extracting data')

        # delete temporary folder
        # shutil.rmtree(temp_folder)

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
    slans_de = SLANSDataExtraction()

    # get simulation data from multiple folders
    dir_path = r'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Data'
    for path in os.listdir(dir_path):
        if os.path.exists(fr"{dir_path}/{path}/shape_space.xlsx"):
            print(path)
            # load shape space
            shape_space = pd.read_excel(fr"{dir_path}/{path}/shape_space.xlsx", "Sheet1", index_col=0)
            shape_space = shape_space.rename('Cavity{}'.format)
            shape_space = shape_space.to_dict(orient='index')
            # ic(shape_space)
            # break
            slans_data_folder = fr"{dir_path}/{path}/SLANS"
            temp_folder = fr"{dir_path}/{path}/temp"
            save_excel = fr"{dir_path}/{path}/results_slans"
            request = '0D Data'
            slans_de.multiple_folders_data_parallel(shape_space, slans_data_folder, 40, 2, 'mm',
                                                    request, save_excel, temp_folder)
            # break
