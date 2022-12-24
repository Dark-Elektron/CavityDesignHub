import json
import os
import subprocess
import sys
import numpy as np

class SLANSPar:
    def __init__(self):
        os.chdir('../../../simulation_codes')
        print(os.getcwd())


        mid_cell_key = input("Please enter mid cell shape key-> ")
        self.f_shift = input("Enter frequency shift (MHz)-> ")
        # mid_cell_key = 46

        # get input dictionary
        self.mid_cell_par_dict = json.load(open('Extracted node_editor\population.json', 'r'))["{}".format(mid_cell_key)]
        self.end_cell_par_dict = json.load(open('Population_C{}\population.json'.format(mid_cell_key), 'r'))

        # get list of input
        self.key_list = []
        population = json.load(open('Population_C{}\population.json'.format(mid_cell_key), 'r'))
        ff_dict = json.load(open('Population_C{}\\ff_population.json'.format(mid_cell_key), 'r'))
        mts = json.load(open('Population_C{}\\mt_score.json'.format(mid_cell_key), 'r'))

        for key, point in population.items():
            if ff_dict[key][0] < 0.98:
                continue
            else:
                # extra code for tm012
                # if np.max(mts[key]['mts']) < 3:
                self.key_list.append(key)
                # else:
                #     continue

        print(len(self.key_list))
        # print(self.mid_cell_par_dict)
        #
        # print(self.end_cell_par_dict)


    def run(self):
        # run slans code
        print(os.getcwd())
        print(self.key_list)

        if len(self.key_list) < 12:
            self.proc_count = len(self.key_list)
        else:
            self.proc_count = 12

        command = ["mpiexec", "-np", "{}".format(self.proc_count), "python", "SLANS/slans_mpi.py",
                   "{}".format(self.key_list), "{}".format(self.mid_cell_par_dict), "{}".format(self.end_cell_par_dict), "{}".format(self.f_shift)]
        # print(command)
        sp = subprocess.run(command, stdout=sys.stdout, stderr=subprocess.STDOUT)
        print("Done")


if __name__ == '__main__':
    slans_run = SLANSPar()
    slans_run.run()
