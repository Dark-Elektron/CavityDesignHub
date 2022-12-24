import json
import os
import time

from SLANS_code.slans_tune import Tune
from SLANS_code.tune_py_Req import TunePy

tune = TunePy()


def init_population(AB_space, ab_space, Ri_space, L, freq):
    population = {}
    freq_dict = {}
    dip_dist = {}
    key = 0
    last_saved_key = -1

    dict_name = '15_48_population.json'
    freq_dict_name = '15_48_freq_dict.json'

    # check for already processed shapes
    response = 'n'
    if os.path.exists(r'C:\Users\sosoho\Dropbox\2D_Codes\ABCI_software\Python_ABCI\15_48_population.json'):
        response = input("Do you want to continue simulation from last key [Y/n]? ")
        if response.strip().lower() == 'y' or response.strip() == '':
            last_saved_key, population, freq_dict = get_latest_key(dict_name, freq_dict_name)
            print("INITIAL POPULATION:: ", population)
    
    for m, AB in enumerate(AB_space):
        for n, ab in enumerate(ab_space):
            for o, Ri in enumerate(Ri_space):
                # edit to check for key later
                if key > last_saved_key:
    
                    # Optimization Constraint:: Req > Ri + B + b
                    # if Req < (Ri + AB + ab):
                    #     Req = Ri + AB + ab + 5
                    Req = Ri + AB + ab + 10

                    # Optimization Constraint, L > A + a
                    if (AB + ab) > L:
                        print("\tCONSTRAINT:: Value set skipped -> A + a > L")
                        continue

                    # Tune cell to get Req.
                    Req, freq = get_req(AB, AB, ab, ab, Req, Ri, L, key)

                    # # Optimization constraint:: alpha > 90 (No reentrant cavity shapes allowed.
                    # if alpha < 90:
                    #     print("\tCONSTRAINT:: Alpha < 90::Values::", AB, AB, ab, ab, Ri, Req)
                    #     continue

                    # # write dictionary dipole
                    population[key] = [AB, AB, ab, ab, Req, Ri, L]
                    freq_dict[key] = freq

                    # update dictionaries
                    print("SLANS_REQ_TUNE:: Writing parameters to dict -> ", population)
                    with open('15_48_population.json', 'w') as file:
                        file.write(json.dumps(population, indent=4, separators=(',', ': ')))
                    with open('15_48_freq_dict.json', 'w') as file:
                        file.write(json.dumps(freq_dict, indent=4, separators=(',', ': ')))

                    print("SLANS_REQ_TUNE:: Done writing")

                key += 1
    
    print("Population::", population)
    # save_to_file()
    return population, dip_dist


def get_req(A, B, a, b, Req, Ri, L, key):

    # B = A/AB
    # b = a/ab
    par_in = [A, B, a, b, Req, Ri, L]
    print("From get Req", par_in)
    Req_l = 170
    Req_r = 170 + 35
    R_eq, freq = tune.window_search(par_in, par_in, Req_l, Req_r, key)

    return R_eq, freq


def get_latest_key(pop_dict_name, freq_dict_name):
    # check if the file exists
    path = r'C:\Users\sosoho\Dropbox\2D_Codes\ABCI_software\Python_ABCI\{}'.format(pop_dict_name)
    freq_path = r'C:\Users\sosoho\Dropbox\2D_Codes\ABCI_software\Python_ABCI\{}'.format(freq_dict_name)

    if os.path.exists(path):

        # load json file
        old_dict = json.load(open(path, 'r'))
        old_freq_dict = json.load(open(freq_path, 'r'))

        # get the latest key
        key_list = [int(f) for f in old_dict.keys()]
        latest_key = max(key_list)
        print("Latest key in dictionary is :: ", latest_key)

        return latest_key, old_dict, old_freq_dict
    else:
        return -1, {}, {}


if __name__ == '__main__':
    AB_space = [14 + 2*x for x in range(18)]
    ab_space = [4.5 + 1.5*x for x in range(7)]
    Ri_space = [50 + 5*x for x in range(7)]

    # AB_space = [48]
    # ab_space = [13.5]
    # Ri_space = [80]

    print(AB_space, ab_space, Ri_space)

    print(len(AB_space)*len(ab_space)*len(Ri_space))

    start = time.time()
    init_population(AB_space, ab_space, Ri_space, 93.5, 801.58)
    stop = time.time()

    print("Run time = ", (stop-start)*1e-3)