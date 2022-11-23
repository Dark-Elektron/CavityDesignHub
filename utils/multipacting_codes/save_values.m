% Program save_values.m
% -------------------------------------------------------------------------
% Save the input arguments to files fieldparam, param, counter_flevels.mat
% and gene_initials.mat.
%
% -------------------------------------------------------------------------
% CALLS TO : check_geodata.m
% 17/04/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% 27/09/00 :                 - emin and emax added to the vector param
% -------------------------------------------------------------------------

% parameters for the field solver
my0    = 4*pi*1e-7; 
eps0   = 8.85418782e-12; 
lambda = 1/(freq*sqrt(eps0*my0));
eta0   = 376.7303134111465;

if gtype == 1                     % constant for streching
  strech = 0;
elseif gtype == 2
  strech = lambda/2;
else  
  strech = lambda/2;              % this may need modification (?)
end

fieldparam = zeros(8,1);
fieldparam(1) = gtype;
save -ascii fieldparam fieldparam

ok = check_geodata;
Z  = 0;
if ok > 0 & gtype > 1
  if gtype == 2
    load geodata.n
    gd = geodata;
  elseif gtype == 3
    load geodatar.n
    gd = geodatar;
  end
  n    = length(gd(:,1));
  zend = max(gd(4:n,2));
  ind  = find(abs(gd(4:n,2)-zend)<=1e-9);
  if length(ind) < 2
    cl = error_message('The ends of the geometry must be vertical.');
  else
    a = min(gd(ind+3,1));
    b = max(gd(ind+3,1));
    if a < 1e-6
      cl = error_message('Can not calculate Z, a = 0.');
    else
      Z = eta0 * log(b/a)/(2*pi);          % impedance of a coaxial line
    end
  end
end

fieldparam = [gtype,freq,epsr,d1/1000,real(R),imag(R),strech,Z]';
save -ascii fieldparam fieldparam

% field levels
flevel = (flmin:flstep:flmax)'*1e3;
if gtype > 1
  flevel = pow2volt(flevel,R,Z);         % change the powers to voltages
end

save counter_flevels flevel

% parameters for the MP analysis
m     = 9.1093879e-31;
q     = 1.6021773e-19;
c     = 2.99792458e8;
param = zeros(7,1);
V0    = 1;                               % intensity of the EM field
%V0    = 2;  error_message('POISTA TÄMÄ - TESTI (save_values.m)');
gamma = c/freq/(2*pi);                   % distance/phase -coefficient
v00   = c*sqrt(1-1/((v0*q/(m*c^2)+1)^2));% initial velocity (relativistic)
%v00   = sqrt(2*v0*q/m);                 % (classical)
ctype = 1;                               % compute counter functions
tol   = 1e-3;                            % tolerance for the ODE solver

param = [freq,V0,gamma,v00,N,ctype,tol,emin,emax]';
save -ascii param param

% grid constant for creating a new grid
ok = check_geodata;
if ok > 0
  if gtype <= 2
    load geodata.n
    geodata(1,1) = d2/1000;
    save -ascii geodata.n geodata
  else
    load geodatal.n
    geodatal(1,1) = d2/1000;
    save -ascii geodatal.n geodatal
  end
end

% parameters for the initial point generator
dx    = dx/1000;                         % dimensions in mm
dc    = dx/10;                           % distance from the nearest corner
alpha = 0;                               % initial velocity angle
dt    = dphi/freq/360;                   % time step
initials = [-dx,dc,alpha,dt,zmin,zmax];
save gene_initials initials
% -----------------------------------------------------------------------

