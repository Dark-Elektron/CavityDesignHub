import json
import os
import time

import numpy as np
import pandas as pd


class FileReader:
    def __init__(self):
        pass

    def to_pd_dataframe(self):
        print("it works")
        #### TAKES IN PYTHON DICTIONARY AS ARGUMENT AND RETURNS PANDAS DATAFRAME
        pass

    def txt_reader(self, filename):
        lines_list = []
        #### OPEN FILE AND UPDATE LINE LIST
        with open(filename, 'r') as f:
            for l in f.readlines():
                line = []
                l = l.strip()
                a = l.split(' ')
                for b in a:
                    if b == '':
                        continue
                    line.append(float(b))
                lines_list.append(line)

        #### CHECK FOR MISSING VALUES IN LIST
        mx_len = 0
        for val in lines_list:
            # GET MAXIMUM LIST LENGTH
            if len(val) > mx_len:
                mx_len = len(val)

        list_update = []
        for val in lines_list:
            if len(val) < mx_len:
                for i in range(0, mx_len - len(val)):
                    val.append(None)
            list_update.append(val)

        #### CONVERT LIST TO PYTHON DICTIONARY
        data_dict = {}
        for key in range(0, len(line)):
            data_dict[key] = [val[key] for val in list_update]

        #### CONVERT DICTIONARY TO PANDAS DATAFRAME
        df = pd.DataFrame.from_dict(data_dict, orient='index').transpose()
        return df

    def excel_reader(self, filename):
        file = pd.ExcelFile(filename)
        # print(file.sheet_names)
        dd = {}
        for sheet_name in file.sheet_names:
            dd[sheet_name] = file.parse(sheet_name)

        return dd

    def svl_reader(self, filename):
        dict = {
              'TITLE': [],
              'CAVITY RADIUS': [],
              'FREQUENCY': [],
              'LENGTH OF WAVE': [],
              'WAVE VALUE': [],
              'QUALITY FACTOR': [],
              'STORED ENERGY': [],
              'TRANSIT TIME': [],
              'EFFECTIVE IMPEDANCE': [],
              'SHUNT IMPEDANCE': [],
              'MAXIMUM MAG. FIELD': [],
              'MAXIMUM ELEC. FIELD': [],
              'ACCELERATION': [],
              'ACCELERATION RATE': [],
              'AVERAGE E.FIELD ON AXIS': [],
              'KM (Emax/Accel.rate)': [],
              'KH (Hmax*Z0/Accel.rate)': [],
        }
        with open(filename, 'r') as f:
            data = f.readlines()
            for i, line in enumerate(data):
                if '*SLANS*' in line:
                    dict['TITLE'].append(line)

                if 'CAVITY RADIUS' in line:
                    dict['CAVITY RADIUS'].append(self._process_line(line))

                if 'FREQUENCY' in line:
                    dict['FREQUENCY'].append(self._process_line(line))

                if 'LENGTH OF WAVE' in line:
                    dict['LENGTH OF WAVE'].append(self._process_line(line))

                if 'WAVE VALUE' in line:
                    dict['WAVE VALUE'].append(self._process_line(line))

                if 'QUALITY FACTOR' in line:
                    dict['QUALITY FACTOR'].append(self._process_line(line))

                if 'STORED ENERGY' in line:
                    dict['STORED ENERGY'].append(self._process_line(line))

                if 'TRANSIT TIME' in line:
                    dict['TRANSIT TIME'].append(self._process_line(line))

                if 'EFFECTIVE IMPEDANCE' in line:
                    dict['EFFECTIVE IMPEDANCE'].append(self._process_line(line))

                if 'SHUNT IMPEDANCE' in line:
                    dict['SHUNT IMPEDANCE'].append(self._process_line(line))

                if 'MAXIMUM MAG. FIELD' in line:
                    dict['MAXIMUM MAG. FIELD'].append(self._process_line(line))

                if 'MAXIMUM ELEC.FIELD' in line:
                    dict['MAXIMUM ELEC. FIELD'].append(self._process_line(line))

                if 'ACCELERATION' in line and not 'RATE' in line:
                    dict['ACCELERATION'].append(self._process_line(line))

                if 'ACCELERATION RATE' in line:
                    dict['ACCELERATION RATE'].append(self._process_line(line))

                if 'AVERAGE E.FIELD ON AXIS' in line:
                    dict['AVERAGE E.FIELD ON AXIS'].append(self._process_line(line))

                if 'KM (Emax/Accel.rate)' in line:
                    dict['KM (Emax/Accel.rate)'].append(self._process_line(line))

                if 'KH (Hmax*Z0/Accel.rate)' in line:
                    dict['KH (Hmax*Z0/Accel.rate)'].append(self._process_line(line))

        return dict

    def top_reader(self):
        pass

    def json_reader(self, dir, header=None):
        df = pd.read_json(dir)

        if header:
            # check if length of header list is same as column length
            if len(header) == len(list(df.columns)):
                df.columns = header
            else:
                print(f'Expected header length of {len(list(df.columns))}, got {len(header)}.')

        return df

    def pam_reader(self):
        pass

    def _process_line(self, line):
        line = line.strip().split(' ')
        res = 0
        for val in line:
            try:
                res = float(val)
                break
            except:
                continue
        return res

    def _combineDict(self, args):

        d1 = json.load(open('Results/population.json', 'r'))
        d2 = json.load(open('Results/population2.json', 'r'))
        d3 = json.load(open('Results/population3.json', 'r'))

        d1.update(d2)
        d1.update(d3)

        with open('population.json', 'w') as file:
            file.write(json.dumps(d1, indent=4, separators=(',', ': ')))

        f1 = json.load(open('Results/freq_diff.json', 'r'))
        f2 = json.load(open('Results/freq_diff2.json', 'r'))
        f3 = json.load(open('Results/freq_diff3.json', 'r'))

        f1.update(f2)
        f1.update(f3)
        with open('freq_diff.json', 'w') as file:
            file.write(json.dumps(f1, indent=4, separators=(',', ': ')))

    def _updateDict(self, d1, d2):
        for key, val_set in d1.items():
            d1[key] += (d2[key])

        return d1


if __name__ == '__main__':
    fr = FileReader()
    # txt = fr.txt_reader(r"D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Data/SLANS/Cavity382/cavity_mm.pam")
    # fr.txt_reader(r"D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Data/SLANS/Cavity1/cavity_mm_7.af")
    # df = fr.excel_reader(r'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Dataframe_full_data.xlsx')
    # print([[x] for x in df["Midcell Data"].iloc[0][1:8]])
    # print(len(df["Midcell Data"]))
    # print([[x] for x in df["Midcell Data"].iloc[0][8:20]])
    # print([[x] for x in df["Midcell Data"].iloc[0][10:11]])

    # # dataframe from json
    # dir = r'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Cavity Population/pseudo_shape_space.json'
    # df = fr.json_reader(dir, ['A', 'B', 'a', 'b', 'Ri', 'L'])
    # for value in df.to_dict(orient='list').values():
    #     print(value)

    # svl reader
    dir = r'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Data\SLANS\Cavity_process_0/cavity_mm.svl'
    d = fr.svl_reader(dir)
    print(d['FREQUENCY'])
