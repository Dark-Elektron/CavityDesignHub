import json
import os
import subprocess
from threading import Thread
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QFileDialog


def run_slans_exe(path, filename=None):
    path = os.path.join(os.getcwd(), path)
    print(path)
    cwd = "{}\SLANS_data".format(os.getcwd())
    if filename:
        t = Thread(target=subprocess.call, args=([path, '{}'.format(filename)],), kwargs={"cwd": cwd})
    else:
        t = Thread(target=subprocess.call, args=(path,), kwargs={"cwd": cwd})

    t.start()

def run_abci_exe(path):
    path = os.path.join(os.getcwd(), path)
    t = Thread(target=subprocess.call, args=(path,))
    t.start()

def show_hide_widget(cb, widget_list):
    if cb.checkState() == 2:
        for widget in widget_list:
            widget.hide()
    else:
        for widget in widget_list:
            widget.show()


def toggle_show_hide(widget):
    if widget.isVisible():
        widget.hide()
    else:
        widget.show()

def treeIndex(widget):
    widget.setCurrentIndex(3)

def get_fid():
    f = QFileDialog()
    foldername = f.getExistingDirectory()
    if foldername:
        fid = list(foldername)[-1]
        return fid

def plot_control(widget):
    # if widget.value() >= 0:
    #     return '{}'.format(widget.value())
    # elif widget.value() == -1:
    #     return '_0'
    if widget.text():
        return widget.text().upper()
    else:
        get_fid()

def get_point_fid(data_point):
    if isinstance(data_point, list):
        pass
    else:
        # cast to float list
        data_point = data_point.split(', ')
        data_point = [float(x) for x in data_point]
        print(data_point)

    # load population
    dict = json.load(open("Extracted node_editor/population.json", "r"))

    if len(data_point) == 3:
        A = data_point[0]
        a = data_point[1]
        Ri = data_point[2]

        for key, v_set in dict.items():
            if A == v_set[0] and a == v_set[2] and Ri == v_set[5]:
                return key
            else:
                continue

    elif len(data_point) > 3:
        A = data_point[0]
        B = data_point[1]
        a = data_point[2]
        b = data_point[3]
        Ri = data_point[4]

        for key, v_set in dict.items():
            if A == v_set[0] and B == v_set[1] and a == v_set[2] and b == v_set[3] and Ri == v_set[5]:
                return key
            else:
                continue

    return 0

def load_result_set(filename, plot_object):
    data_set = json.load(open(filename, 'r'))

    i = 0
    print(data_set)
    for fid, data in data_set.items():
        if i == 0:
            plot_object.graph_results(fid)
        else:
            plot_object.add_plots(fid)

        i += 1

def read_parameters(ui):
    print("Got here")
    data_set = json.load(open("Cavity Population/population.json", 'r'))

    fid = ui.le_Parameters_Fid.text().upper().strip()
    print(fid)
    try:
        fid = "{}".format(int(fid))
        data = data_set[fid]

        ui.dsb_Aeq_M.setValue(data[0])
        ui.dsb_Beq_M.setValue(data[1])
        ui.dsb_ai_M.setValue(data[2])
        ui.dsb_bi_M.setValue(data[3])
        ui.dsb_Req_M.setValue(data[4])
        ui.dsb_Ri_M.setValue(data[5])
        ui.dsb_L_M.setValue(data[6])

        A, B, a, b, Req, Ri, L = data[0], data[1], data[2], data[3], data[4], data[5], data[6]

        return [A, B, a, b, Req, Ri, L]

    except:
        try:
            path = os.path.join(os.getcwd(), "SLANS_data/Cavity{}/cst_parameters_mid.txt".format(fid))

            var_dict = {"Req_M": 0,
                        "ri_M": 0,
                        'L_M': 0,
                        'Aeq_M': 0,
                        'Beq_M': 0,
                        'ai_M': 0,
                        'bi_M': 0}
            key_list = ["Req_M", "ri_M", 'L_M', 'Aeq_M', 'Beq_M', 'ai_M', 'bi_M']

            with open(path, 'r') as f:
                data = f.readlines()
                for i, line in enumerate(data):
                    var_dict[key_list[i]] = process_line(line)

            Req, Ri, L, A, B, a, b = list(var_dict.values())

            ui.dsb_Aeq_M.setValue(A)
            ui.dsb_Beq_M.setValue(B)
            ui.dsb_ai_M.setValue(a)
            ui.dsb_bi_M.setValue(b)
            ui.dsb_Req_M.setValue(Req)
            ui.dsb_Ri_M.setValue(Ri)
            ui.dsb_L_M.setValue(L)

            return [A, B, a, b, Req, Ri, L]

        except:
            print("node_editor not yet found")
            return [0, 0, 0, 0, 0, 0, 0]


def process_line(line):
    line = line.strip().split('=')
    try:
        result = float(line[-1])
        return result

    except ValueError:
        print("Check file! Cannot get the correct variables.")

def updateDict(dict_dir):
    dict = json.load(open(dict_dir, 'r'))
    population = json.load(open("Extracted node_editor/population.json", 'r'))
    for key, val_set in dict.items():

        population[key] += (dict[key])

    return population

def conv_dict_dataframe(dict, trim):
    data = np.array(list(dict.values()))
    
    data = data[:, trim]
    print("it gets here2")
    dataframe = pd.DataFrame(data=data, index=[i for i in range(np.shape(data)[0])], columns=["A", "a", "Req", "Ri"])
    return dataframe

def plotSNS(win):
    ui = win.ui
    # get dict
    dict = updateDict(get_filename())
    trim = ui.le_Trim_Dataframe.text()
    if trim == "":
        trim = (0, 2, 4, 5)

    print("it gets here1")
    # convert dict to dataframe
    dataframe = conv_dict_dataframe(dict, trim)

    g = sns.pairplot(dataframe, hue='Ri', corner=True, palette="bright", diag_kind="hist", kind='scatter')
    print("it gets here3")
    print(type(win.plt.fig))
    g.fig = win.plt.fig
    print("it gets here5")
    g.axes = win.plt.fig.add_subplot(211)#win.plt.axes
    print("it gets here6")
    win.plt.draw()
    # plt.show()

def get_filename():
    options = QFileDialog.Options()
    # options |= QFileDialog.DontUseNativeDialog
    filename, _ = QFileDialog.getOpenFileName(None, "Select .json", "",
                                                  "All Files (*);;Text Files (*.txt)", options=options)

    print(filename)
    if filename != "":
        if ".json" in filename:
            print(filename)
            print("it gets here4")
            return filename

def pause_thread(t):
    t.set()

def resume_thread(t):
    t.clear()

def stop_thread(t):
    t.stop()

