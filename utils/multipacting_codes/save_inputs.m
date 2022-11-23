% Program save_inputs.m
% -----------------------------------------------------------------------
% Save and check the validity of the input parameters.
%
% -----------------------------------------------------------------------
% CALLS TO : save_values.m, generate_initials.m, error_message.m, 
%            check_frequency.m, check_epsilon.m, check_gridconstant1.m,
%            check_refcoeff.m, check_gridconstant2.m, check_numimpacts.m,
%            check_inivelocity.m, chech_phasestep.m, check_spacestep.m,
%            check_spacez.m, check_flevelstep.m, check_flevels.m
% 06/03/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% 19/05/00 :                 - check for inputs added
% 27/09/00 :                 - inputs emin and emax added
% -----------------------------------------------------------------------

function save_inputs

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

% save the inputs 
save_values;

% and generate initial points
ok = generate_initials(1);

load fieldparam
if fieldparam(1) < 3
  load counter_initials
  nci = length(initials(:,1));
else  
  load counter_initialsr
  load counter_initialsl
  ncir = length(initialsr);
  ncil = length(initialsl);
end  

load counter_flevels
nfl = length(flevel);

cl = 0;
if ok > 0
  cl = error_message(['  ' num2str(nfl) ' field levels.']);
  if fieldparam(1) < 3
    cl=error_message(['  ' num2str(nci) ' initial points.']);
  else
    cl=error_message(['  ' num2str(ncir) ' initial points on the cold side.']);
    cl=error_message(['  ' num2str(ncil) ' initial points on the warm side.']);
  end
  cl=error_message('Input arguments saved.');
end

if cl == 1
  clear_window;
end  
% -----------------------------------------------------------------------
