import json


class UQ:
    def __init__(self):
        pass

    def uq(key, shape, qois, n_cells, n_modules, n_modes, f_shift, bc, parentDir, projectDir):
        err = False
        result_dict_slans = {}
        slans_obj_list = qois
        for o in qois:
            result_dict_slans[o] = {'expe': [], 'stdDev': []}

        # EXAMPLE: p_true = np.array([1, 2, 3, 4, 5]).T
        p_true = shape['IC'][0:5]
        # ic(p_true)
        rdim = len(p_true)  # How many variabels will be considered as random in our case 5
        degree = 1

        #  for 1D opti you can use stroud5 (please test your code for stroud3 less quadrature nodes 2rdim)
        flag_stroud = 1
        if flag_stroud == 1:
            nodes, weights, bpoly = quad_stroud3(rdim, degree)
            nodes = 2. * nodes - 1.
        elif flag_stroud == 2:
            nodes, weights, bpoly = quad_stroud3(rdim, degree)  # change to stroud 5 later
            nodes = 2. * nodes - 1.
        else:
            ic('flag_stroud==1 or flag_stroud==2')

        #  mean value of geometrical parameters
        p_init = np.zeros(np.shape(p_true))

        no_parm, no_sims = np.shape(nodes)
        ic(nodes)
        delta = 0.005  # or 0.1

        Ttab_val_f = []
        print_('3')
        sub_dir = fr'Cavity{key}'  # the simulation runs at the quadrature points are saved to the key of mean value run
        for i in range(no_sims):
            skip = False
            p_init[0] = p_true[0] * (1 + delta * nodes[0, i])
            p_init[1] = p_true[1] * (1 + delta * nodes[1, i])
            p_init[2] = p_true[2] * (1 + delta * nodes[2, i])
            p_init[3] = p_true[3] * (1 + delta * nodes[3, i])
            p_init[4] = p_true[4] * (1 + delta * nodes[4, i])

            par_mid = np.append(p_init, shape['IC'][5:]).tolist()
            par_end = par_mid
            print_(par_mid)

            # perform checks on geometry
            ok = perform_geometry_checks(par_mid, par_end)
            print_("OK", ok)
            if not ok:
                err = True
                break
            fid = fr'{key}_Q{i}'

            # check if folder exists and skip if it does
            print_(fr'{projectDir}\SimulationData\SLANS\Cavity{key}\Cavity{fid}')
            if os.path.exists(fr'{projectDir}\SimulationData\SLANS\Cavity{key}\Cavity{fid}'):
                skip = True
                # ic("Skipped: ", fid, fr'{projectDir}\SimulationData\ABCI\Cavity{key}\Cavity{fid}')

            # skip analysis if folder already exists.
            if not skip:
                #  run model using SLANS or CST
                # # create folders for all keys
                slans_geom.createFolder(fid, projectDir, subdir=sub_dir)
                try:
                    slans_geom.cavity(n_cells, n_modules, par_mid, par_end, par_end,
                                      n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                      parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)
                except KeyError:
                    slans_geom.cavity(n_cells, n_modules, par_mid, par_end, par_end,
                                      n_modes=n_modes, fid=fid, f_shift=f_shift, bc=bc, beampipes=shape['BP'],
                                      parentDir=parentDir, projectDir=projectDir, subdir=sub_dir)

            filename = fr'{projectDir}\SimulationData\SLANS\Cavity{key}\Cavity{fid}\cavity_{bc}.svl'
            if os.path.exists(filename):
                params = fr.svl_reader(filename)
                norm_length = 2 * n_cells * shape['IC'][5]
                ic(n_cells, norm_length)
                qois_result = get_qoi_value(params, slans_obj_list, n_cells, norm_length)
                print_(qois_result)
                # sometimes some degenerate shapes are still generated and the solver returns zero
                # for the objective functions, such shapes are considered invalid
                for objr in qois_result:
                    if objr == 0:
                        # skip key
                        err = True
                        break

                tab_val_f = qois_result

                Ttab_val_f.append(tab_val_f)
            else:
                err = True

        # # add original point
        # filename = fr'{projectDir}\SimulationData\SLANS\Cavity{key}\cavity_33.svl'
        # params = fr.svl_reader(filename)
        # obj_result, tune_result = get_objectives_value(params, slans_obj_list)
        # tab_val_f = obj_result
        # Ttab_val_f.append(tab_val_f)

        print_("Error: ", err)
        ic(Ttab_val_f)
        # import matplotlib.pyplot as plt
        if not err:
            v_expe_fobj, v_stdDev_fobj = weighted_mean_obj(np.atleast_2d(Ttab_val_f), weights)
            # ic(v_expe_fobj, v_stdDev_fobj)
            # append results to dict
            ic(v_expe_fobj, v_stdDev_fobj)
            for i, o in enumerate(slans_obj_list):
                result_dict_slans[o]['expe'].append(v_expe_fobj[i])
                result_dict_slans[o]['stdDev'].append(v_stdDev_fobj[i])

                # pdf = normal_dist(np.sort(np.array(Ttab_val_f).T[i]), v_expe_fobj[i], v_stdDev_fobj[i])
                # plt.plot(np.sort(np.array(Ttab_val_f).T[i]), pdf)

            # plt.show()

            with open(fr"{projectDir}\SimulationData\SLANS\Cavity{key}\uq.json", 'w') as file:
                file.write(json.dumps(result_dict_slans, indent=4, separators=(',', ': ')))
        else:
            print_(fr"There was a problem running UQ analysis for Cavity{key}")

    def get_qoi_value(d, obj, n_cells, norm_length):
        Req = d['CAVITY RADIUS'][n_cells - 1] * 10  # convert to mm
        Freq = d['FREQUENCY'][n_cells - 1]
        E_stored = d['STORED ENERGY'][n_cells - 1]
        Rsh = d['SHUNT IMPEDANCE'][n_cells - 1]  # MOhm
        Q = d['QUALITY FACTOR'][n_cells - 1]
        Epk = d['MAXIMUM ELEC. FIELD'][n_cells - 1]  # MV/m
        Hpk = d['MAXIMUM MAG. FIELD'][n_cells - 1]  # A/m
        # Vacc = dict['ACCELERATION'][0]
        Eavg = d['AVERAGE E.FIELD ON AXIS'][n_cells - 1]  # MV/m
        Rsh_Q = d['EFFECTIVE IMPEDANCE'][n_cells - 1]  # Ohm

        Vacc = np.sqrt(
            2 * Rsh_Q * E_stored * 2 * np.pi * Freq * 1e6) * 1e-6
        # factor of 2, remember circuit and accelerator definition
        # Eacc = Vacc / (374 * 1e-3)  # factor of 2, remember circuit and accelerator definition
        Eacc = Vacc / (norm_length * 1e-3)  # for 1 cell factor of 2, remember circuit and accelerator definition
        Epk_Eacc = Epk / Eacc
        Bpk_Eacc = (Hpk * 4 * np.pi * 1e-7) * 1e3 / Eacc

        d = {
            "Req": Req,
            "freq": Freq,
            "Q": Q,
            "E": E_stored,
            "R/Q": 2 * Rsh_Q,
            "Epk/Eacc": Epk_Eacc,
            "Bpk/Eacc": Bpk_Eacc
        }

        objective = []

        # append objective functions
        for o in obj:
            if o in d.keys():
                objective.append(d[o])

        return objective
