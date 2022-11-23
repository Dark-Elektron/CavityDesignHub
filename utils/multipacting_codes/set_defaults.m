% Program set_defaults.m
% -----------------------------------------------------------------------
% Set the defaults values for the input arguments.
%
% -----------------------------------------------------------------------
% CALLS TO : save_values.m, error_message.m
% 06/03/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% 27/09/00 :                 - parameters emin and emax added
% -----------------------------------------------------------------------

function set_defaults

figure(2)

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
set(fobj,'String',10);
d1   = 10;

% reflection coeff
fobj = findobj('Tag','RefCoeffRe');
set(fobj,'String',-1);

fobj = findobj('Tag','RefCoeffIm');
set(fobj,'String',0);
R    = -1;

% grid constant for the MP analysis
fobj = findobj('Tag','GridCons2');
set(fobj,'String',2.5);
d2   = 2.5;

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
set(fobj,'String',-100);
zmin = -100;

% space step for initial sites
fobj = findobj('Tag','SpaceStep');
set(fobj,'String',3);
dx   = 3;

% max z for initial sites
fobj = findobj('Tag','Maxz');
set(fobj,'String',100);
zmax = 100;

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
set(fobj,'String',40e3);
flmax = 40e3;

save_values;

error_message('Default values for the inputs are set.');
% -----------------------------------------------------------------------

