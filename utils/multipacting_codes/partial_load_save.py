# Clear all variables
import numpy as np
clear('all')
# Create dummy variables to test with
a = randi(10,1000,1)
b = randi(10,1000,1)
c = randi(10,1000,1)
# Save variables a & b into a v4 type matfile
save('test_v4.mat','-v4')
# Clear variables in workspace
clear('all')
# Get list of variables in "test_v4.mat"
vars_in_testv4 = whos('-file','multipac/H.mat')
# Loop over variables names and load into workspace one by one
for i in np.arange(1,len(vars_in_testv4)+1).reshape(-1):
    # Load variable into workspace from "test_v4.mat"
    scipy.io.loadmat('multipac/H.mat',vars_in_testv4(i).name)

# Save variables into a new matfile with v7 format
save('test_v7.mat',vars_in_testv4.name,'-v7')