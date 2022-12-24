import ast
from sys import argv
from mpi4py import MPI
from termcolor import colored
from analysis_modules.eigenmode.SLANS.slans_geom_par import SLANSGeometry

slans_geom = SLANSGeometry()

file_color = 'cyan'

def print_(*arg):
    print(colored(f'\t\t\t{arg}', file_color))

def run_mpi():
    # print_("\tRUNNING MPI")
    # get variables from argv

    n_cells = ast.literal_eval(u'{}'.format(argv[1]))
    n_modules = ast.literal_eval(u'{}'.format(argv[2]))
    shape_space = ast.literal_eval(u'{}'.format(argv[3]))
    f_shift = ast.literal_eval(u'{}'.format(argv[4]))
    n_modes = ast.literal_eval(u'{}'.format(argv[5]))
    parentDir = ast.literal_eval(u'{}'.format(argv[6]))
    projectDir = ast.literal_eval(u'{}'.format(argv[7]))

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

    # print_("Process {} got {}".format(rank, proc_list))

    for key in proc_keys_list:# get parameters from key
        try:
            # run slans code
            slans_geom.cavity(n_cells, n_modules, shape_space[key]['MC'], shape_space[key]['LC'], shape_space[key]['RC'],
                              n_modes=n_modes, fid="{}".format(key), f_shift=f_shift, beampipes=shape_space[key]['BP'],
                              parentDir=parentDir, projectDir=projectDir)

        except Exception as e:
            print_("EXCEPTION CAUGHT:: ", e)
            print_("\t\tPROCESS {}, KEY {} COULD NOT COMPLETE THE CODE".format(rank, key))
            pass

if __name__ == '__main__':
    # print("\tRUNNING MPI from main")
    run_mpi()