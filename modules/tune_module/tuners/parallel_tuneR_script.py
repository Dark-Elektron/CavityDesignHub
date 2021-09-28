import ast
import os
import stat
# from subprocess import call
from math import floor
from sys import argv

# import mpi4py
# mpi4py.rc.recv_mprobe = False
from mpi4py import MPI
from utils.file_reader import FileReader
fr = FileReader()
import numpy as np
from termcolor import colored
from SLANS_code.slans_geom_par import SLANSGeometry
import shutil
from distutils.dir_util import copy_tree
from genetic_algorithm import Evolution
import time as t

evol = Evolution()
slans_geom = SLANSGeometry()

file_color = 'cyan'

DEBUG = False
# DEBUG = True

def print_(*arg):
    if DEBUG: print(colored(fr'\t{arg}', file_color, on_color='on_white'))


def run_mpi():
    print_("\t\tRUNNING MPI")
    # get variables from argv
    # pseudo_shape_space = ast.literal_eval(f'{argv[1]}')
    pseudo_shape_space_name = argv[1]
    resume = argv[2]
    last_key = argv[3]
    parentDir = argv[4]

    # get dictionary from json file
    dirc = fr'{parentDir}/{pseudo_shape_space_name}'
    pseudo_shape_space = fr.json_reader(dirc)

    # update/resume dictionary: get only dictionary of cavities not yet processed
    trimmed_dict = {}
    for key, val in pseudo_shape_space.items():
        if int(key) >= int(last_key):
            trimmed_dict[key] = val

    pseudo_shape_space = trimmed_dict
    # print(pseudo_shape_space)

    # print(pseudo_shape_space)
    # print(type(pseudo_shape_space))

    # print_(f'mid: {mid_cell_par}, ss: {pseudo_shape_space}, resume: {resume}, parentDir: {parentDir}')

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    print_(f'{rank}, {size}')

    # copy needed files with only main process to avoid data race. Barrier is placed so that other proccesses wait
    if rank == 0:
        print_("Rank copying necessary files")
        for i in range(size):
            # remove folders to avoid reading of previous data
            overwriteFolder(i, parentDir)
            copyFiles(i, parentDir)
        print_("done copying files")

    comm.barrier()

    # overwriteFolder(rank, parentDir)
    # copyFiles(rank, parentDir)

    # share pseudo_shape_space amongst processors
    keys_list = list(pseudo_shape_space.keys())
    pseudo_shape_space_len = len(keys_list)
    chunk_size = floor(pseudo_shape_space_len/size)

    # print_(f'ss_len: {pseudo_shape_space_len}, chunk_size: {chunk_size}')

    if rank < size - 1:
        print_(rank)
        proc_keys = keys_list[rank*chunk_size: chunk_size*(rank + 1)]
        pseudo_shape_space_proc = {f'{key}': pseudo_shape_space[key] for key in proc_keys}
    else:
        print_(rank)
        proc_keys = keys_list[rank*chunk_size:]
        pseudo_shape_space_proc = {f'{key}': pseudo_shape_space[key] for key in proc_keys}

    # print_(f'ss_proc {rank}: {pseudo_shape_space_proc}')

    # try:
    # run slans code
    evol.tune(pseudo_shape_space_proc,,  #, last_key=last_key This would have to be tested again #val2
    # except Exception as e:
    #     print_(f"PROCESS {rank} COULD NOT COMPLETE THE CODE -> {e}")

def overwriteFolder(invar, parentDir):
    print_("IT's in overwrite")
    path = f"{parentDir}\Data\SLANS\Cavity_process_{invar}"
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)

    # if os.path.exists(path):
    #     # print("\t\t\tremoved directory:: ", path)
    #     cond = True
    #     n = 0
    #     while cond:
    #         # print("Attempt:: ", n)
    #         try:
    #             shutil.rmtree(path, onerror=remove_readonly)
    #             # print("\t\tdone removing")
    #             os.mkdir(path)
    #             cond = False
    #         except Exception as e:
    #             print_(f'Exception encountered trying to overwrite files -> {e}')
    #             n += 1
    #             continue
    # else:
    #     print_("Creating new folder, first run!")
    #     os.mkdir(path)


def copyFiles(invar, parentDir):
    src = f"{parentDir}\SLANS_exe"
    dst = f"{parentDir}\Data\SLANS\Cavity_process_{invar}\SLANS_exe"
    # print(src)
    # print(dst)
    # call(['cp', '-a', src, dst])
    copy_tree(src, dst)


def remove_readonly(func, path, excinfo):
    os.chmod(path, stat.S_IWRITE)
    func(path)


if __name__ == '__main__':
    run_mpi()