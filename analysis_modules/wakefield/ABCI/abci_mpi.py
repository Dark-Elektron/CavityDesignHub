import ast
# from subprocess import call
import time
from sys import argv
from mpi4py import MPI
# from ABCI.abci_geom_par import ABCIGeometry
from analysis_modules.wakefield.ABCI.abci_geometry import ABCIGeometry
from termcolor import colored

abci_geom = ABCIGeometry()

file_color = 'cyan'

def print_(*arg):
    print(colored(f'\t\t\t{arg}', file_color))

def run_mpi():
    print_("\tRUNNING MPI")
    # get variables from argv

    n_cells = ast.literal_eval(u'{}'.format(argv[1]))
    n_modules = ast.literal_eval(u'{}'.format(argv[2]))
    shape_space = ast.literal_eval(u'{}'.format(argv[3]))
    wakefield_parameters = ast.literal_eval(u'{}'.format(argv[4]))
    mesh_parameters = ast.literal_eval(u'{}'.format(argv[5]))
    parentDir = ast.literal_eval(u'{}'.format(argv[6]))
    projectDir = ast.literal_eval(u'{}'.format(argv[7]))

    MROT, MT, NFS, UBT, bunch_length = wakefield_parameters
    DDR_SIG, DDZ_SIG = mesh_parameters

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    # share process by rank
    shape_space_len = len(shape_space)
    share = round(shape_space_len/size)

    keys = list(shape_space.keys())

    if rank < size:
        proc_keys_list = keys[rank*share:rank*share + share]
    else:
        proc_keys_list = keys[rank*share:]

    print_("Process {} got {}".format(rank, proc_keys_list))

    for key in proc_keys_list:
        # get parameters from key
        try:
            # run slans code
            # print_(f'Key: {key}')
            start = time.time()
            abci_geom.cavity(n_cells, n_modules, shape_space[key]['MC'], shape_space[key]['LC'], shape_space[key]['RC'],
                                  fid=key, MROT=MROT, MT=MT, NFS=NFS, UBT=UBT, bunch_length=bunch_length,
                                  DDR_SIG=DDR_SIG, DDZ_SIG=DDZ_SIG, parentDir=parentDir, projectDir=projectDir)
            end = time.time()
            print_(f"Proc {rank} is done running {key}, {end-start}s")

        except Exception as e:
            print_("EXCEPTION CAUGHT:: ", e)
            print_("\t\tPROCESS {}, KEY {} COULD NOT COMPLETE THE CODE".format(rank, key))
            pass

if __name__ == '__main__':
    # print_("\tRUNNING MPI from main")
    run_mpi()