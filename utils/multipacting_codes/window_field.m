% Function program window_field(s)
% ------------------------------------------------------------------------
% Runs the field solver for computing the EM fields in a window geometry.
%
% ------------------------------------------------------------------------
% CALLS TO : error_message.m, make_model_window.m, plot_mesh.m, 
%            find_resonance.m, calculate_fields.m, plot_FEM_fields.m,
% 13/04/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% -------------------------------------------------------------------------

function window_field(s)

if nargin < 1
  s = 1;
end

load fieldparam
freq    = fieldparam(2);
releps  = fieldparam(3);
gridcons= fieldparam(4);

if releps == 1
  error_message('Relative permittivity == 1?');
end

%load geodatal.n
gridcons1 = geodatal(1,1);
  
% generate the mesh
error_message('Generating the mesh.');

make_model_window(gridcons,releps,s);

% plot the mesh
figure(3)
plot_mesh(s);

error_message('Searching the resonances. Blue area for streching.');
[k1,k2] = find_resonance(freq);

mu0   = 4*pi*1e-7; 
eps0  = 8.85418782e-12; 
freq1 = k1/(2*pi*sqrt(mu0*eps0));
freq2 = k2/(2*pi*sqrt(mu0*eps0));
err1  = abs(freq1-freq)/freq*100
err2  = abs(freq2-freq)/freq*100
if max(err1,err2) > 1
  error_message('Error in eigen frequency more than 1%.');
end  

% compute, plot and save the fields
calculate_fields(1,gridcons1,s);

figure(4)
plot_FEM_fields(0,gridcons1,s);
% ---------------------------------------------------------------------

