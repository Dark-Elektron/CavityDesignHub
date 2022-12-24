import ast
import multiprocessing as mp
import time
from sys import argv
from utils.file_reader import FileReader
# from simulation_codes.ABCI.abci_geom_par import ABCIGeometry
from modules.wakefield.ABCI.abci_geometry import ABCIGeometry
from termcolor import colored

abci_geom = ABCIGeometry()
fr = FileReader()

file_color = 'cyan'


def print_(*arg):
    print(colored(f'\t\t\t{arg}', file_color))


def run_MP():
    print_("\tRUNNING MPI")
    # get variables from argv

    n_cells = ast.literal_eval(u'{}'.format(argv[1]))
    n_modules = ast.literal_eval(u'{}'.format(argv[2]))
    shape_space_name = f'{argv[3]}'
    wakefield_parameters = ast.literal_eval(u'{}'.format(argv[4]))
    mesh_parameters = ast.literal_eval(u'{}'.format(argv[5]))
    proc_count = int(argv[6])
    parentDir = f'{argv[7]}'
    projectDir = f'{argv[8]}'

    MROT, MT, NFS, UBT, bunch_length = wakefield_parameters
    DDR_SIG, DDZ_SIG = mesh_parameters

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
            if p < proc_count - 1:
                proc_keys_list = keys[p * share:p * share + share]
            else:
                proc_keys_list = keys[p * share:]

            processor_shape_space = {}
            for key, val in shape_space.items():
                if proc_keys_list[0] <= int(key) <= proc_keys_list[-1]:
                    processor_shape_space[key] = val
            # print(f'Processor {p}: {processor_shape_space}')

            service = mp.Process(target=run_sequential, args=(n_cells, n_modules, processor_shape_space,
                                                              MROT, MT, NFS, UBT, bunch_length,
                                                              DDR_SIG, DDZ_SIG,
                                                              parentDir, projectDir
                                                              ))
            service.start()
            processes.append(service)
            # print("Done")

        except Exception as e:
            print_("Exception in run_MP::", e)


def run_sequential(n_cells, n_modules, processor_shape_space,
                   MROT=0, MT=4, NFS=10000, UBT=50, bunch_length=20,
                   DDR_SIG=0.1, DDZ_SIG=0.1,
                   parentDir=None, projectDir=None):
    for key, shape in processor_shape_space.items():
        try:
            # # create folders for all keys
            # slans_geom.createFolder(key)

            # run slans code
            start_time = time.time()
            try:
                abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC_R'],
                                 fid=key, MROT=MROT, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                 DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir)
            except:
                abci_geom.cavity(n_cells, n_modules, shape['IC'], shape['OC'], shape['OC'],
                                 fid=key, MROT=MROT, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                 DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir)

            print_(f'Cavity {key}. Time: {time.time() - start_time}')
        except Exception as e:
            print(f'Error in abci_mpi_mp:: run_sequential -> {e}')


if __name__ == '__main__':
    # print_("\tRUNNING MPI from main")
    run_MP()
