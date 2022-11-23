% Clear all variables
clear all;

% Create dummy variables to test with
a = randi(10,1000,1);
b = randi(10,1000,1);
c = randi(10,1000,1);

% Save variables a & b into a v4 type matfile
save('test_v4.mat', '-v4');

% Clear variables in workspace
clear all;

% Get list of variables in "test_v4.mat"
vars_in_testv4 = whos('-file', 'multipac/H.mat');

% Loop over variables names and load into workspace one by one
for i = 1: length(vars_in_testv4)
    
    % Load variable into workspace from "test_v4.mat"
    load('multipac/H.mat', vars_in_testv4(i).name);
end

% Save variables into a new matfile with v7 format
save("test_v7.mat", vars_in_testv4.name, '-v7');




