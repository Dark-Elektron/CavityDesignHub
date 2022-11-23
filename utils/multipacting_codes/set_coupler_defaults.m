% Program set_coupler_defaults.m
% -----------------------------------------------------------------------
% Set the defaults values for the input arguments due to a coupler.
% (A segment of a coaxial line).
% -----------------------------------------------------------------------
% CALLS TO : save_values.m, 
% 11/05/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% 27/09/00 :                 - parameters emin and emax added
% -----------------------------------------------------------------------

function set_coupler_defaults

figure(2)

% type of the geometry
fobj  = findobj('Tag','GeoType');
set(fobj,'Value',2);
gtype = 2;

% frequency
fobj = findobj('Tag','Frequency');
set(fobj,'String',1.3);
freq = 1.3e9;

% relative epsilon of the window
fobj = findobj('Tag','Epsilon');
set(fobj,'String',1);
epsr = 1;

% grid constant for the field solver
fobj = findobj('Tag','GridCons1');
set(fobj,'String',4);
d1   = 4;

% reflection coeff
fobj = findobj('Tag','RefCoeffRe');
set(fobj,'String',-1);

fobj = findobj('Tag','RefCoeffIm');
set(fobj,'String',0);
R    = -1;

% grid constant for computing and plotting the fields
fobj = findobj('Tag','GridCons2');
set(fobj,'String',2);
d2   = 2;

% number of impacts
fobj = findobj('Tag','NumOfIm');
set(fobj,'String',20);
N    = 20;

% initial velocity
fobj = findobj('Tag','Velocity');
set(fobj,'String',2);
v0   = 2;

% min impact energy
fobj = findobj('Tag','Emmin');
set(fobj,'String',0);
emin  = 0;

% max impact energy
fobj = findobj('Tag','Emmax');
set(fobj,'String',1e6);
emax  = 1e6;

% phase step for initial sites
fobj = findobj('Tag','PhaseStep');
set(fobj,'String',5);
dphi = 5;

% min z for initial sites
fobj = findobj('Tag','Minz');
set(fobj,'String',0.055);
zmin = 0.055;

% space step for initial sites
fobj = findobj('Tag','SpaceStep');
set(fobj,'String',1);
dx   = 1;

% max z for initial sites
fobj = findobj('Tag','Maxz');
set(fobj,'String',0.06);
zmax = 0.06;

% min field level
fobj  = findobj('Tag','MinFl');
set(fobj,'String',50);
flmin = 50;

% field level step
fobj  = findobj('Tag','StepFl');
set(fobj,'String',5);
flstep = 5;

% max field level
fobj  = findobj('Tag','MaxFl');
set(fobj,'String',1000);
flmax = 1000;

save_values;
% -----------------------------------------------------------------------

