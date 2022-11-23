% Program save_parameters.m
% -------------------------------------------------------------------------
% Save the input arguments to file mpgui_inputs.mat.
%
% -------------------------------------------------------------------------
% CALLS TO : check_geodata.m
% 15/12/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% -------------------------------------------------------------------------

function save_parameters

figure(2)

% type of the geometry
fobj  = findobj('Tag','GeoType');
gtype = get(fobj,'Value');

% frequency
fobj = findobj('Tag','Frequency');
freq = str2num(get(fobj,'String')) * 1e9;
ok   = check_frequency(freq); if ok == 0, return; end

% relative epsilon of the window
fobj = findobj('Tag','Epsilon');
epsr = str2num(get(fobj,'String'));
ok   = check_epsilon(epsr,gtype); if ok == 0, return; end

% grid constant for the field solver
fobj = findobj('Tag','GridCons1');
d1   = str2num(get(fobj,'String'));
ok   = check_gridconstant1(d1); if ok == 0, return; end

% reflection coeff
fobj = findobj('Tag','RefCoeffRe');
Rre  = str2num(get(fobj,'String'));
ok   = check_refcoeff1(Rre); if ok == 0, return; end

fobj = findobj('Tag','RefCoeffIm');
Rim  = str2num(get(fobj,'String'));
ok   = check_refcoeff1(Rim); if ok == 0, return; end

R    = Rre + i*Rim;
ok   = check_refcoeff(R,gtype); if ok == 0, return; end

% grid constant for the MP analysis
fobj = findobj('Tag','GridCons2');
d2   = str2num(get(fobj,'String'));
ok   = check_gridconstant2(d2); if ok == 0, return; end

% number of impacts
fobj = findobj('Tag','NumOfIm');
N    = str2num(get(fobj,'String'));
ok   = check_numimpacts(N); if ok == 0, return; end

% initial velocity
fobj = findobj('Tag','Velocity');
v0   = str2num(get(fobj,'String'));
ok   = check_inivelocity(v0); if ok == 0, return; end

% minimum impact energy
fobj = findobj('Tag','Emmin');
emin = str2num(get(fobj,'String'));

% maximum impact energy
fobj = findobj('Tag','Emmax');
emax = str2num(get(fobj,'String'));
ok   = check_energy(emin,emax); if ok == 0, return; end

% phase step for initial sites
fobj = findobj('Tag','PhaseStep');
dphi = str2num(get(fobj,'String'));
ok   = check_phasestep(dphi); if ok == 0, return; end

% z min for initial sites
fobj = findobj('Tag','Minz');
zmin = str2num(get(fobj,'String'));

% space step for initial sites
fobj = findobj('Tag','SpaceStep');
dx   = str2num(get(fobj,'String'));
ok   = check_spacestep(dx); if ok == 0, return; end

% max z for initial sites
fobj = findobj('Tag','Maxz');
zmax = str2num(get(fobj,'String'));
ok   = check_spacez(zmin,zmax,gtype); if ok == 0, return; end

% min field level
fobj  = findobj('Tag','MinFl');
flmin = str2num(get(fobj,'String'));

% field level step
fobj   = findobj('Tag','StepFl');
flstep = str2num(get(fobj,'String'));
ok     = check_flevelstep(flstep); if ok == 0, return; end

% max field level
fobj  = findobj('Tag','MaxFl');
flmax = str2num(get(fobj,'String'));
ok    = check_flevels(flmin,flmax); if ok == 0, return; end

% save the inputs parameters to the file mpgui_inputs.mat
save mpgui_inputs gtype freq epsr d1 R d2 N v0 emin emax dphi zmin dx ...
     zmax flmin flstep flmax

error_message('Parameters of the Input Window saved.');
% -----------------------------------------------------------------------
