import os
import subprocess


class Dakota:
    def __init__(self, folder, name):
        self.folder = folder
        self.name = name

    def write_input_file(self, kind, descriptors, responses, **kwargs):
        # check if folder exists, if not, create folder
        if not os.path.exists(os.path.join(self.folder, self.name)):
            try:
                os.mkdir(os.path.join(self.folder, self.name))
            except:
                print("Could not create folder. Make sure target location exists.")

        with open(os.path.join(self.folder, self.name, f'{self.name}.in'), 'w') as f:
            self.environment(f)
            self.method(f)
            self.variables(f, kind, descriptors, **kwargs)
            self.interface(f, responses)
            self.responses(f, responses)

    def environment(self, f):
        f.write('environment\n')
        f.write('\ttabular_data\n')
        f.write("\t\ttabular_data_file = 'sim_result_table.dat'\n")
        f.write('\tresults_output\n')
        f.write("\t\tresults_output_file = 'result_output_file.dat'\n")

        f.write('\n')

    def method(self, f, method='polynomial_chaos', **kwargs):
        f.write('method\n')
        f.write(f'\t{method}\n')
        f.write("\t\texport_expansion_file ='expansion_file.dat'\n")
        f.write("\t\tcubature_integrand = 3\n")
        f.write("\t\tsamples_on_emulator = 10000\n")
        f.write("\t\tseed = 12347\n")
        f.write("\t\tvariance_based_decomp interaction_order = 1\n")

        f.write('\n')

    def variables(self, f, kind, descriptors, **kwargs):
        """

        Parameters
        ----------
        f: File
            File
        kind: {'uniform_uncertain', 'normal_uncertain', 'beta_uncertain'}
            Type of distribution of the variables
        descriptors: list, ndarray
            Uncertain variable names
        kwargs: kwargs
            Other relevant arguments eg. {means: [], 'std_deviations': [], 'lower_bounds': [], 'upper_bounds': []}


        Returns
        -------

        """
        f.write("variables\n")
        f.write(f"\t{kind} = {len(descriptors)}\n")
        f.write("\tdescriptors       =   " + '\t\t\t'.join(['"'+descriptor+'"' for descriptor in descriptors]) + '\n')

        if 'means' in kwargs.keys():
            assert len(descriptors) == len(kwargs['means'])
            f.write("\tmeans      =   " + '\t\t\t'.join([str(mean) for mean in kwargs['means']]) + '\n')

        if 'std_deviations' in kwargs.keys():
            assert len(descriptors) == len(kwargs['std_deviations'])
            f.write("\tstd_deviations      =   " + '\t\t\t'.join([str(std) for std in kwargs['std_deviations']]) + '\n')

        if 'lower_bounds' in kwargs.keys():
            assert len(descriptors) == len(kwargs['lower_bounds'])
            f.write("\tlower_bounds      =   " + '\t\t\t'.join([str(lb) for lb in kwargs['lower_bounds']]) + '\n')

        if 'upper_bounds' in kwargs.keys():
            assert len(descriptors) == len(kwargs['upper_bounds'])
            f.write("\tupper_bounds      =   " + '\t\t\t'.join([str(ub) for ub in kwargs['upper_bounds']]) + '\n')

        f.write('\n')

    def interface(self, f, responses, num_procs=1):
        f.write("interface\n")
        f.write("#\tcommon options\n")
        f.write("#\tfork\n")
        f.write("\tparameters_file = 'params.in'\n")
        f.write("\tresults_file    = 'results.out'\n")
        f.write(f"\tsystem asynchronous evaluation_concurrency = {num_procs}\n")
        f.write(f"\tanalysis_driver = 'python pycall.py {len(responses)} True'\n")
        f.write("#\tparameters_file = 'params.in'\n")
        f.write("#\tresults_file    = 'results.out'\n")
        f.write("#\tfile_tag\n")
        f.write("#\tfile_save\n")
        f.write("#\taprepro\n")

        f.write('\n')

    def responses(self, f, responses):
        f.write("responses\n")
        f.write(f"\tresponse_functions = {len(responses)}\n")
        f.write("\tno_gradients\n")
        f.write("\tno_hessians\n")

        f.write('\n')

    def run_analysis(self, nodes_only=False):
        cwd = os.path.join(self.folder, self.name)
        subprocess.run(['dakota', '-i', f'{os.path.join(self.folder, self.name, f"{self.name}.in")}', '-o', f'{os.path.join(self.folder, self.name, f"{self.name}.out")}'],
                       cwd=cwd, shell=True)


if __name__ == '__main__':
    folder = r'D:\Dropbox\CavityDesignHub\analysis_modules\uq\dakota_scripts'
    name = 'test_dakotapy'
    dakota = Dakota(folder, name)
    kind = 'uniform_uncertain'
    parameters = ['A', 'B']
    resp = ['f1', 'f2']
    lower_bounds = [0, 0]
    upper_bounds = [1, 2]
    dakota.write_input_file(kind, descriptors=parameters, responses=resp, lower_bounds=lower_bounds, upper_bounds=upper_bounds)

    dakota.run_analysis()