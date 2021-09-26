import ast
import multiprocessing as mp
import time
from sys import argv

import numpy as np
from termcolor import colored
from simulation_codes.SLANS.slans_geom_par import SLANSGeometry
from utils.file_reader import FileReader

# evol = Evolution()
slans_geom = SLANSGeometry()
fr = FileReader()


file_color = 'cyan'

def print_(*arg):
    print(colored(f'\t\t\t{arg}', file_color))

def run_MP():
    print_("\t\tRUNNING MP")

    # get variables from argv
    n_cells = ast.literal_eval(u'{}'.format(argv[1]))
    n_modules = ast.literal_eval(u'{}'.format(argv[2]))
    shape_space_name = f'{argv[3]}'
    f_shift = ast.literal_eval(u'{}'.format(argv[4]))
    n_modes = ast.literal_eval(u'{}'.format(argv[5]))
    proc_count = int(argv[6])
    parentDir = argv[7]
    projectDir = argv[8]

    # get dictionary from json file
    dirc = fr'{projectDir}\Cavities\{shape_space_name}'
    shape_space = fr.json_reader(dirc)

    # split shape_space for different processes/ MPI share process by rank
    keys = list(shape_space.keys())
    shape_space_len = len(keys)
    share = round(shape_space_len / proc_count)

    processes = []
    for p in range(proc_count):
        try:
            if p < proc_count-1:
                proc_keys_list = keys[p * share:p * share + share]
            else:
                proc_keys_list = keys[p * share:]

            processor_shape_space = {}
            for key, val in shape_space.items():
                if proc_keys_list[0] <= int(key) <= proc_keys_list[-1]:
                    processor_shape_space[key] = val
            # print(f'Processor {p}: {processor_shape_space}')

            service = mp.Process(target=run_sequential, args=(n_cells, n_modules, processor_shape_space, n_modes, f_shift, parentDir, projectDir))
            service.start()
            processes.append(service)
            # print("Done")

        except Exception as e:
            print_("Exception in run_MP::", e)

def run_sequential(n_cells, n_modules, processor_shape_space, n_modes, f_shift, parentDir, projectDir):
    for key, shape in processor_shape_space.items():
        try:
            # # create folders for all keys
            # slans_geom.createFolder(key)

            # run slans code
            start_time = time.time()
            try:
                slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                  n_modes=n_modes, fid=f"{key}", f_shift=f_shift, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir)
            except:
                slans_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                  n_modes=n_modes, fid=f"{key}", f_shift=f_shift, beampipes=shape['BP'],
                                  parentDir=parentDir, projectDir=projectDir)

            print_(f'Done with Cavity {key}. Time: {time.time() - start_time}')
        except Exception as e:
            print(f'Error in slans_mpi_mp:: run_sequential -> {e}')


if __name__ == '__main__':
    # print("\tRUNNING MPI from main")
    run_MP()