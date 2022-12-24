import os

import matplotlib.pyplot as plt
from utils.file_reader import FileReader
fr = FileReader()

class Monitor:
    def generate_data(self):
        # check if convergence json already exists
        if os.path.exists(''):
            self.update_data()
        else:
            pass

    def update_data(self):
        data = fr.json_reader('')

        # get list of folders not empty
        non_empty_folders_key = []

        # compare list with json key
        for k in non_empty_folders_key:
            if k not in data.keys():
                # process the folder

                # update json
                data[k] = []
                # with open()


    def get_convergence_data(self, monitor_list):
        min = {'k_loss_long': 2, 'k_loss_trans': 2, 'df_M0D1': 2}
        data = fr.excel_reader(fr'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\modules\data_module\combined_data.xlsx')
        data = data['First_Batch']
        data_m = {}
        data_conv = {}
        # get only data to be monitored
        for entry in monitor_list:
            data_m[entry] = [x / data[entry][0] for x in data[entry]]
            data_conv[entry] = []

        # convert to convergence plot
        for par, vals in data_m.items():
            for val in vals:
                if val < min[par]:
                    data_conv[par].append(val)
                    min[par] = val
                else:
                    data_conv[par].append(min[par])

        return data_conv
