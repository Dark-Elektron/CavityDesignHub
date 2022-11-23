% Function program show_initials.m
% --------------------------------------------------------------------
% Generates and plots the initial sites for the MP analysis.
%
% --------------------------------------------------------------------
% CALLS TO : check_param.m, generate_initials.m, plot_initials.m
%            error_message.m
% 07/03/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% ---------------------------------------------------------------------

function show_initials

figure(2)

% frequency
fobj = findobj('Tag','Frequency');
freq = str2num(get(fobj,'String')) * 1e9;

% get the new values for the parameters
% phase step for initial sites
fobj = findobj('Tag','PhaseStep');
dphi = str2num(get(fobj,'String'));

% z min for initial sites
fobj = findobj('Tag','Minz');
zmin = str2num(get(fobj,'String'));

% space step for initial sites
fobj = findobj('Tag','SpaceStep');
dx   = str2num(get(fobj,'String'));

% max z for initial sites
fobj = findobj('Tag','Maxz');
zmax = str2num(get(fobj,'String'));

dx    = dx/1000;                         % dimensions in mm
dc    = dx/10;                           % distance from the nearest corner
alpha = 0;                               % initial velocity angle
dt    = dphi/freq/360;                   % time step
initials = [-dx,dc,alpha,dt,zmin,zmax];

save gene_temp initials

% type of the geometry
fobj  = findobj('Tag','GeoType');
gtype = get(fobj,'Value');

ok = 1;
ok = check_param;
if gtype <= 2
  fileis = exist('geodata.n');
else
  fileis1 = exist('geodatal.n');
  fileis2 = exist('geodataw.n');
  fileis3 = exist('geodatar.n');
  fileis  = fileis1*fileis2*fileis3;
end  
if fileis == 0
  error_message('Geometry file is missing.');
  ok = 0;
end

if ok > 0
  generate_initials(0);
  plot_initials(0)
  error_message('To save the inputs press Save.');
  error_message('                              ');
end  
% ---------------------------------------------------------------------

