import ast
import multiprocessing as mp
import os
import stat
# from subprocess import call
from sys import argv

from termcolor import colored
from analysis_modules.eigenmode.SLANS.slans_geom_par import SLANSGeometry
import shutil
from distutils.dir_util import copy_tree
from analysis_modules.cavity_analysis.analysis_codes import Analysis
from utils.file_reader import FileReader

analysis = Analysis()
slans_geom = SLANSGeometry()
fr = FileReader()

file_color = 'cyan'

# DEBUG = True
DEBUG = True

def print_(*arg):
    if DEBUG: print(colored(fr'\t{arg}', file_color, on_color='on_white'))


def run_MP():
    print_("\t\tRUNNING MP")

    # get variables from argv
    proc_count = int(argv[1])
    pseudo_shape_space_name = f'{argv[2]}'
    resume = argv[3]
    last_key = argv[4]
    parentDir = argv[5]
    projectDir = argv[6]
    tuner = argv[7]
    tune_variable = argv[8]
    iter_set = ast.literal_eval(argv[9])
    cell_type = argv[10]

    # get dictionary from json file
    dirc = fr'{pseudo_shape_space_name}'
    pseudo_shape_space = fr.json_reader(dirc)

    # update/resume dictionary: get only dictionary of cavities not yet processed
    trimmed_dict = {}
    for key, val in pseudo_shape_space.items():
        if int(key) >= int(last_key):
            trimmed_dict[key] = val

    pseudo_shape_space = trimmed_dict

    # split shape_space for different processes/ MPI share process by rank
    keys = list(pseudo_shape_space.keys())
    shape_space_len = len(keys)
    share = round(shape_space_len / proc_count)

    processes = []
    for p in range(proc_count):
        try:
            if p < proc_count:
                proc_keys_list = keys[p * share:p * share + share]
            else:
                proc_keys_list = keys[p * share:]

            overwriteFolder(p, projectDir)
            copyFiles(p, parentDir, projectDir)

            processor_shape_space = {}
            for key, val in pseudo_shape_space.items():
                if proc_keys_list[0] <= int(key) <= proc_keys_list[-1]:
                    processor_shape_space[key] = val
            print(f'Processor {p}: {processor_shape_space}')

            service = mp.Process(target=run_sequential, args=(processor_shape_space, resume, p, parentDir, projectDir, tuner,
                                                              tune_variable, iter_set, cell_type))
            service.start()
            processes.append(service)

        except Exception as e:
            print_("Exception in run_MP::", e)

def run_sequential(pseudo_shape_space_proc, resume, p, parentDir, projectDir, tuner, tune_variable, iter_set, cell_type):
    analysis.GSSEC(pseudo_shape_space_proc, parentDir, projectDir, resume=resume, proc=p, tuner=tuner,
                   tune_variable=tune_variable, iter_set=iter_set, cell_type=cell_type) #, last_key=last_key This would have to be tested again #val2

def overwriteFolder(invar, projectDir):
    print_("IT's in overwrite")
    path = f"{projectDir}\SimulationData\SLANS\Cavity_process_{invar}"
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)

def copyFiles(invar, parentDir, projectDir):
    src = fr"{parentDir}\exe\SLANS_exe"
    dst = fr"{projectDir}\SimulationData\SLANS\Cavity_process_{invar}\SLANS_exe"
    copy_tree(src, dst)

def remove_readonly(func, path, excinfo):
    os.chmod(path, stat.S_IWRITE)
    func(path)

if __name__ == '__main__':
    run_MP()