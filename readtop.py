import time

filename = r'D:\Dropbox\2D_Codes\ABCI_software\Python_ABCI\Data\ABCI\Cavity_3\Cavity_MROT_0.top'
frame_objects = {}
frame_titles_objects = {}
frame_count = 0
new_line_count = 0

t = time.time()
plot_decorations = [r'Cavity Shape Input', r'Cavity Shape Used', r'Wake Potentials', r'Real Part of Longitudinal Impedance', r'Imaginary Part of Longitudinal Impedance',
                    'Frequency Spectrum of Loss Factor', r'Loss Factor Spectrum Integrated up to F', r'Real Part of Long. + Log Impedance', r'Imaginary Part of Long. + Log Impedance',
                    r'Spectrum of Long. + Log Loss Factor', r'Long. + Log Factor Integrated up to F',  r'Real Part of Azimuthal Impedance',
                    r'Imaginary Part of Azimuthal Impedance', r'Real Part of Transverse Impedance', r'Imaginary Part of Transverse Impedance']

with open(filename, 'r') as f:
    for line in f:
        if 'NEW FRAME' in line:
            # add frame to frame_objects
            if frame_count > 0:
                frame_titles_objects[frame_count-1] = frame_title

                # select suitable title
                for decor in plot_decorations:
                    if decor in frame_title[0]+frame_title[1]:
                        frame_objects[decor] = line_objects
                        break

            # reset lists
            line_objects = []
            frame_title = []
            new_frame = True


            frame_count += 1

        if new_frame:
            new_frame = False

        # add titles to frame_title
        if 'TITLE' in line or 'MORE' in line:
            frame_title.append(line)

        if 'SET LIMITS' in line:
            x, y = [], []

        if 'JOIN' in line:
            new_line_object = True
        else:
            new_line_object = False

        if new_line_object:
            # add x and y to line objects
            line_objects.append([x, y])

            # reset x and y
            x, y = [], []

        if len(line.strip().split()) == 2 and line.strip().split()[0] != 'NEW' and line.strip().split()[0] != 'JOIN':
            new_line_object = False
            ll = [float(a) for a in line.strip().split()]
            x.append(ll[0])
            y.append(ll[1])

    # EOF append last list
    frame_titles_objects[frame_count - 1] = frame_title
    # select suitable title
    for decor in plot_decorations:
        if decor in frame_title[0]+frame_title[1]:
            frame_objects[decor] = line_objects
            break

    new = {}
    for key, lines in frame_objects.items():
        x, y = [], []
        for line in lines:
            if len(line[0]) > 2:
                x += line[0]
                y += line[1]

        new[key] = [x, y]

    print('runtime: ', time.time() - t)

    # import matplotlib.pyplot as plt
    # xy = new[list(new.keys())[2]]
    # plt.plot(xy[0], xy[1])
    # plt.show()




