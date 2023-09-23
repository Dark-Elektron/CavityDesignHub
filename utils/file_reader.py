import json
import pandas as pd


class FileReader:
    def __init__(self):
        pass

    def to_pd_dataframe(self):
        print("it works")
        # TAKES IN PYTHON DICTIONARY AS ARGUMENT AND RETURNS PANDAS DATAFRAME
        pass

    def txt_reader(self, filename, delimeter):
        lines_list = []
        # OPEN FILE AND UPDATE LINE LIST
        with open(filename, 'r') as f:
            for l in f.readlines():
                line = []
                l = l.strip()
                a = l.split(delimeter)
                for b in a:
                    if b == '':
                        continue
                    line.append(float(b))
                lines_list.append(line)

        # CHECK FOR MISSING VALUES IN LIST
        mx_len = 0
        for val in lines_list:
            # GET MAXIMUM LIST LENGTH
            if len(val) > mx_len:
                mx_len = len(val)

        list_update = []
        for val in lines_list:
            if len(val) < mx_len:
                for i in range(0, mx_len - len(val)):
                    val.append(None)
            list_update.append(val)

        # CONVERT LIST TO PYTHON DICTIONARY
        data_dict = {}
        for key in range(0, len(line)):
            data_dict[f"{key}"] = [val[key] for val in list_update]

        # CONVERT DICTIONARY TO PANDAS DATAFRAME
        df = pd.DataFrame.from_dict(data_dict, orient='index').transpose()
        return df

    def excel_reader(self, filename):
        file = pd.ExcelFile(filename)
        # print(file.sheet_names)
        dd = {}
        for sheet_name in file.sheet_names:
            dd[sheet_name] = file.parse(sheet_name)

        return dd

    def svl_reader(self, filename):
        dict = {
            'TITLE': [],
            'CAVITY RADIUS': [],
            'LENGTH': [],
            'FREQUENCY': [],
            'LENGTH OF WAVE': [],
            'WAVE VALUE': [],
            'QUALITY FACTOR': [],
            'STORED ENERGY': [],
            'TRANSIT TIME': [],
            'EFFECTIVE IMPEDANCE': [],
            'SHUNT IMPEDANCE': [],
            'MAXIMUM MAG. FIELD': [],
            'MAXIMUM ELEC. FIELD': [],
            'ACCELERATION': [],
            'ACCELERATION RATE': [],
            'AVERAGE E.FIELD ON AXIS': [],
            'KM (Emax/Accel.rate)': [],
            'KH (Hmax*Z0/Accel.rate)': [],
        }
        with open(filename, 'r') as f:
            data = f.readlines()
            for i, line in enumerate(data):
                if '*SLANS*' in line:
                    dict['TITLE'].append(line)

                if 'CAVITY RADIUS' in line:
                    dict['CAVITY RADIUS'].append(self._process_line(line, 'CAVITY RADIUS'))
                    dict['LENGTH'].append(self._process_line(line, 'LENGTH'))

                if 'FREQUENCY' in line:
                    dict['FREQUENCY'].append(self._process_line(line, 'FREQUENCY'))

                if 'LENGTH OF WAVE' in line:
                    dict['LENGTH OF WAVE'].append(self._process_line(line, 'LENGTH OF WAVE'))

                if 'WAVE VALUE' in line:
                    dict['WAVE VALUE'].append(self._process_line(line, 'WAVE VALUE'))

                if 'QUALITY FACTOR' in line:
                    dict['QUALITY FACTOR'].append(self._process_line(line, 'QUALITY FACTOR'))

                if 'STORED ENERGY' in line:
                    dict['STORED ENERGY'].append(self._process_line(line, 'STORED ENERGY'))

                if 'TRANSIT TIME' in line:
                    dict['TRANSIT TIME'].append(self._process_line(line, 'TRANSIT TIME'))

                if 'EFFECTIVE IMPEDANCE' in line:
                    dict['EFFECTIVE IMPEDANCE'].append(self._process_line(line, 'EFFECTIVE IMPEDANCE'))

                if 'SHUNT IMPEDANCE' in line:
                    dict['SHUNT IMPEDANCE'].append(self._process_line(line, 'SHUNT IMPEDANCE'))

                if 'MAXIMUM MAG. FIELD' in line:
                    dict['MAXIMUM MAG. FIELD'].append(self._process_line(line, 'MAXIMUM MAG. FIELD'))

                if 'MAXIMUM ELEC.FIELD' in line:
                    dict['MAXIMUM ELEC. FIELD'].append(self._process_line(line, 'MAXIMUM ELEC.FIELD'))

                if 'ACCELERATION' in line and not 'RATE' in line:
                    dict['ACCELERATION'].append(self._process_line(line, 'ACCELERATION'))

                if 'ACCELERATION RATE' in line:
                    dict['ACCELERATION RATE'].append(self._process_line(line, 'ACCELERATION RATE'))

                if 'AVERAGE E.FIELD ON AXIS' in line:
                    dict['AVERAGE E.FIELD ON AXIS'].append(self._process_line(line, 'AVERAGE E.FIELD ON AXIS'))

                if 'KM (Emax/Accel.rate)' in line:
                    dict['KM (Emax/Accel.rate)'].append(self._process_line(line, 'KM (Emax/Accel.rate)'))

                if 'KH (Hmax*Z0/Accel.rate)' in line:
                    dict['KH (Hmax*Z0/Accel.rate)'].append(self._process_line(line, 'KH (Hmax*Z0/Accel.rate)'))

        return dict

    def sv2_reader(self, filename):
        dict = {
            'Number of azimuth variation': [],
            'Frequency': [],
            'Length of wave': [],
            'Wave value': [],
            'Quality factor': [],
            'Stored energy': [],
            'Transverse impedance/Q': [],
            'Eff.transverce impedance': [],
            'Eff.shunt tran.impedance': [],
            'Enmax/Hrmax': [],
            'Enmax': [],
            'Zemax(cm), Remax(cm)': [],
            'Htmax': [],
            'Zhmax(cm), Rhmax(cm)': [],
            'Htmax/Heff': [],
            'Hrmin(A/m), Z(cm)': [],
            'Hrmax(A/m), Z(cm)': []
        }

        with open(filename, 'r') as f:
            data = f.readlines()
            for i, line in enumerate(data):
                if 'Number of azimuth variation' in line:
                    dict['Number of azimuth variation'].append(self._process_line(line, 'Number of azimuth variation'))

                if 'Frequency' in line:
                    dict['Frequency'].append(self._process_line(line, 'Frequency'))

                if 'Length of wave' in line:
                    dict['Length of wave'].append(self._process_line(line, 'Length of wave'))

                if 'Wave value' in line:
                    dict['Wave value'].append(self._process_line(line, 'Wave value'))

                if 'Quality factor' in line:
                    dict['Quality factor'].append(self._process_line(line, 'Quality factor'))

                if 'Stored energy' in line:
                    dict['Stored energy'].append(self._process_line(line, 'Stored energy'))

                if 'Transverse impedance/Q' in line:
                    dict['Transverse impedance/Q'].append(self._process_line(line, 'Transverse impedance/Q'))

                if 'Eff.transverce impedance' in line:
                    dict['Eff.transverce impedance'].append(self._process_line(line, 'Eff.transverce impedance'))

                if 'Eff.shunt tran.impedance' in line:
                    dict['Eff.shunt tran.impedance'].append(self._process_line(line, 'Eff.shunt tran.impedance'))

                if 'Enmax/Hrmax' in line:
                    dict['Enmax/Hrmax'].append(self._process_line(line, 'Enmax/Hrmax'))

                if 'Enmax' in line:
                    dict['Enmax'].append(self._process_line(line, 'Enmax'))

                if 'Zemax(cm), Remax(cm)' in line and not 'RATE' in line:
                    dict['Zemax(cm), Remax(cm)'].append(self._process_line(line, 'Zemax(cm)'))
                    dict['Zemax(cm), Remax(cm)'].append(self._process_line(line, 'Remax(cm)'))

                if 'Htmax' in line:
                    dict['Htmax'].append(self._process_line(line, 'Htmax'))

                if 'Zhmax(cm), Rhmax(cm)' in line:
                    dict['Zhmax(cm), Rhmax(cm)'].append(self._process_line(line, 'Zhmax(cm)'))
                    dict['Zhmax(cm), Rhmax(cm)'].append(self._process_line(line, 'Rhmax(cm)'))

                if 'Htmax/Heff' in line:
                    dict['Htmax/Heff'].append(self._process_line(line, 'Htmax/Heff'))

                if 'Hrmin(A/m), Z(cm)' in line:
                    dict['Hrmin(A/m), Z(cm)'].append(self._process_line(line, 'Hrmin(A/m)'))
                    dict['Hrmin(A/m), Z(cm)'].append(self._process_line(line, 'Z(cm)'))

                if 'Hrmax(A/m), Z(cm)' in line:
                    dict['Hrmax(A/m), Z(cm)'].append(self._process_line(line, 'Hrmax(A/m)'))
                    dict['Hrmax(A/m), Z(cm)'].append(self._process_line(line, 'Z(cm)'))
        return dict

    def top_reader(self):
        pass

    @staticmethod
    def json_reader(filename, header=None):

        df = pd.read_json(filename)

        if header:
            # check if length of header list is same as column length
            if len(header) == len(list(df.columns)):
                df.columns = header
            else:
                print(f'Expected header length of {len(list(df.columns))}, got {len(header)}.')

        try:  # this is to try to return a dataframe of cavity variables if the json file is a shape space
            t = df[:][:1].loc["IC"]
            dd = t.apply(pd.Series)
            dd = dd.set_axis(["A", "B", "ai", "bi", "Ri", "L", "Req", "alpha"], axis=1, inplace=False)

            return dd  # dd is a dataframe whose columns are the cavity geometry variables
        except:
            return df

    def pam_reader(self):
        pass

    def _process_line(self, line, request):
        # select substring from index
        line = line[line.index(request):]
        line = line.strip().split(" ")
        # print(line)
        res = 0
        for val in line:
            try:
                res = float(val)
                break
            except:
                continue
        return res

    def _combineDict(self, args):

        d1 = json.load(open('Results/population.json', 'r'))
        d2 = json.load(open('Results/population2.json', 'r'))
        d3 = json.load(open('Results/population3.json', 'r'))

        d1.update(d2)
        d1.update(d3)

        with open('population.json', 'w') as file:
            file.write(json.dumps(d1, indent=4, separators=(',', ': ')))

        f1 = json.load(open('Results/freq_diff.json', 'r'))
        f2 = json.load(open('Results/freq_diff2.json', 'r'))
        f3 = json.load(open('Results/freq_diff3.json', 'r'))

        f1.update(f2)
        f1.update(f3)
        with open('freq_diff.json', 'w') as file:
            file.write(json.dumps(f1, indent=4, separators=(',', ': ')))

    def _updateDict(self, d1, d2):
        for key, val_set in d1.items():
            d1[key] += (d2[key])

        return d1


if __name__ == '__main__':
    fr = FileReader()

    # txt reader
    durr = r"D:\Dropbox\CavityDesignHub\SampleProject\Cavities\test.json"

    d = fr.json_reader(durr)
    print(type(d))
    print(d)
    # expand df.tags into its own dataframe
    t = d[:][:1].loc["IC"]
    dd = t.apply(pd.Series)
    dd = dd.set_axis(["A", "B", "a", "b", "Ri", "L", "Req", "alpha"], axis=1, inplace=False)
