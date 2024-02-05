import os


class Dakota:
    def __init__(self, folder, name):
        self.folder = folder
        self.name = name

    def write_input_file(self, kind, descriptors, responses):
        # check if folder exists, if not, create folder
        if not os.path.exists(self.folder):
            try:
                os.mkdir(os.path.join(self.folder, self.name))
            except:
                print("Could not create folder. Make sure target location exists.")

        with open(os.path.join(self.folder, self.name, f'{self.name}.in')) as f:
            self.environment(f)
            self.method(f)
            self.variables(f, kind, descriptors)
            self.interface(f)
            self.responses(f, responses)

    def environment(self, f):
        f.write('environment')
        f.write('\ttabular_data')
        f.write("\t\ttabular_data_file = 'sim_result_table.dat'")
        f.write('\tresults_output')
        f.write("\t\tresults_output_file = 'result_output_file.dat'")

    def method(self, f):
        f.write('method')
        f.write('\tpolynomial_chaos')
        f.write("\t\texport_expansion_file ='expansion_file.dat'")
        f.write("\t\tcubature_integrand = 3")
        f.write("\t\tsamples_on_emulator = 10000")
        f.write("\t\tseed = 12347")
        f.write("\t\tvariance_based_decomp interaction_order = 1")

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
        f.write("variables")
        f.write(f"\t{kind} = {len(descriptors)}")
        f.write("\tdescriptors       =   " + '\t'.join([descriptor for descriptor in descriptors]))

        if 'means' in kwargs.keys():
            assert len(descriptors) == len(kwargs['means'])
            f.write("\tmeans      =   " + '\t'.join([lb for lb in kwargs['means']]))

        if 'std_deviations' in kwargs.keys():
            assert len(descriptors) == len(kwargs['std_deviations'])
            f.write("\tstd_deviations      =   " + '\t'.join([lb for lb in kwargs['std_deviations']]))

        if 'lower_bounds' in kwargs.keys():
            assert len(descriptors) == len(kwargs['lower_bounds'])
            f.write("\tlower_bounds      =   " + '\t'.join([lb for lb in kwargs['lower_bounds']]))

        if 'upper_bounds' in kwargs.keys():
            assert len(descriptors) == len(kwargs['upper_bounds'])
            f.write("\tupper_bounds      =   " + '\t'.join([ub for ub in kwargs['upper_bounds']]))

    def interface(self, f, num_procs=1):
        f.write("interface")
        f.write("#\tcommon options")
        f.write("#\tfork,")
        f.write("\tparameters_file = 'params.in'")
        f.write("\tresults_file    = 'results.out'")
        f.write(f"\tsystem asynchronous evaluation_concurrency = {num_procs}")
        f.write("\tanalysis_driver = 'python3 pycall.py'")
        f.write("#\tparameters_file = 'params.in'")
        f.write("#\tresults_file    = 'results.out'	")
        f.write("#\tfile_tag")
        f.write("#\tfile_save")
        f.write("#\taprepro")

    def responses(self, f, responses):
        f.write("responses")
        f.write(f"\tresponse_functions = {len(responses)}")
        f.write("\tno_gradients")
        f.write("\tno_hessians")

if __name__ == '__main__':
    folder = r'D:\Dropbox\CavityDesignHub\analysis_modules\uq\dakota_scripts'
    name = 'test_dakotapy'
    dakota = Dakota(folder, name)
    kind = 'normal_uncertain'
    parameters = ['A', 'B']
    resp = ['f1']
    dakota.write_input_file(kind, descriptors=parameters, responses=resp)