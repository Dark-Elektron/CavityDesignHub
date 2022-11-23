% Program set_cavity_defaults.m
% -----------------------------------------------------------------------
% Set the defaults values for the input arguments due to the TESLA cavity.
%
% -----------------------------------------------------------------------
% CALLS TO : save_values.m
% 11/05/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% -----------------------------------------------------------------------

function set_cavity_defaults

% type of the geometry
fobj  = findobj('Tag','GeoType');
set(fobj,'Value',1);
gtype = 1;

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
set(fobj,'String',5);
d1   = 5;

% reflection coeff
fobj = findobj('Tag','RefCoeffRe');
set(fobj,'String',-1);

fobj = findobj('Tag','RefCoeffIm');
set(fobj,'String',0);
R    = -1;

% grid constant for computing and plotting the fields
fobj = findobj('Tag','GridCons2');
set(fobj,'String',3);
d2   = 3;

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
set(fobj,'String',-0.04);
zmin = -0.04;

% space step for initial sites
fobj = findobj('Tag','SpaceStep');
set(fobj,'String',3);
dx   = 3;

% max z for initial sites
fobj = findobj('Tag','Maxz');
set(fobj,'String',0);
zmax = 0;

% min field level
fobj  = findobj('Tag','MinFl');
set(fobj,'String',1e3);
flmin = 1e3;

% field level step
fobj  = findobj('Tag','StepFl');
set(fobj,'String',0.5e3);
flstep = 0.5e3;

% max field level
fobj  = findobj('Tag','MaxFl');
set(fobj,'String',60e3);
flmax = 60e3;

save_values;
% -----------------------------------------------------------------------
