import json
import os
import subprocess
import sys

DEBUG = True
class ABCIPar:
    def __init__(self):
        os.chdir('..')
        print(os.getcwd())
        self.proc_count = 8
        mid_cell_key = 46

        # get input dictionary
        self.mid_cell_par = json.load(open('Extracted Data\population.json', 'r'))["{}".format(mid_cell_key)]
        self.end_cell_par_dict = json.load(open('Population_C46\population.json', 'r'))

        # get list of input
        self.key_list = []
        population = json.load(open('Population_C46\population.json', 'r'))
        ff_dict = json.load(open('Population_C46\\ff_population.json', 'r'))

        for key, point in population.items():
            if ff_dict[key][0] < 0.98:
                continue
            else:
                self.key_list.append(key)

        if DEBUG: print(len(self.key_list))
        if DEBUG: print(self.mid_cell_par)

        if DEBUG: print(self.end_cell_par_dict)

    def run(self):
        # run abci code
        if DEBUG: print(os.getcwd())
        command = ["mpiexec", "-np", "{}".format(self.proc_count), "python", "ABCI/abci_mpi.py",
                   "{}".format(self.key_list), "{}".format(self.mid_cell_par), "{}".format(self.end_cell_par_dict)]
        if DEBUG: print(command)
        sp = subprocess.run(command, stdout=sys.stdout, stderr=subprocess.STDOUT)
        if DEBUG: print("Done")


if __name__ == '__main__':
    abci_run = ABCIPar()
    abci_run.run()