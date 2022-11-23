% Function program P0 = peak_coupler_field(wall)
% ------------------------------------------------------------------------
% Computes the rf power in a coupler and window geometry. The power is 
% determined at the right wall (on the cavity or cold side).
%
% INPUT wall : 0 - a cavity with electric walls
%              1 - a coupler or window with electric walls
%              2 - a coupler or window with magnetic walls
% OUTPUT  P0 : rf power, [W]
% -----------------------------------------------------------------------
% CALLS TO : field_bin, gaussleg.m
% 09/05/00 : Pasi Ylä-Oijala - Rolf Nevalinna Institute
% -----------------------------------------------------------------------

function P0 = peak_coupler_field(wall)

% load mesh
load mesh
load fieldparam
gtype = fieldparam(1);
Z     = fieldparam(8);
if gtype == 1
  P0 = 0;
else
  P0 = zeros(1,3);
end  

if wall == 0
  !copy /y kama0.mat kama.mat  
  % compute the field at the right end
  load geodata.n
else
  if wall == 1
    !copy /y kama1.mat kama.mat
  else 
    !copy /y kama2.mat kama.mat
  end
  % compute the field at the right end
  if gtype == 2
    load geodata.n
  elseif gtype == 3
    load geodatar.n
    geodata = geodatar;
  end
end

n   = length(geodata(:,1));
zmax= max(geodata(4:n,2));
ind = find(abs(geodata(4:n,2)-zmax)<=1e-9);
if length(ind) < 2
  error_message('The ends of the geometry must be vertical.');
  return;
else
  rmin = min(geodata(ind+3,1));
  rmax = max(geodata(ind+3,1));
  r = rmin:(rmax-rmin)/19:rmax; 
  r = r(:);
  z = zmax*ones(size(r)) - 1e-9;
end
rind = 1:length(r); rind = rind(:);
zind = 1:length(z); zind = zind(:);
alue = 0; % Uses the whole mesh (except those elements marked with 
          %  0:s = the stretch area ) 

save zr z -v4
save zr r -v4 -append
save zr alue -v4 -append
save zr rind -v4 -append
save zr zind -v4 -append

% compute the field at generated points
!Multipac fields

load Er
load Ez
load H

ee   = sqrt(Er.^2+Ez.^2);
eta0 = 376.7303134111465;
hh   = abs(H)*eta0;

[x,w] = gaussleg(min(r),max(r),20);
if wall == 2            % magnetic walls  
  Ex = interp1(r,Er,x,'cubic');
  U  = sum(w.*abs(Ex));
  if sum(Er) >= 0
    if gtype == 1
      P0    = -max(ee);                       % field level      
    else
      P0(1) = -U;                             % voltage
      P0(2) = -U^2/(8*Z);                     % power
      P0(3) = -max(ee);                       % field level
    end
  else
    if gtype == 1
      P0    =  max(ee);                       % field level
    else
      P0(1) =  U;                             % voltage 
      P0(2) =  U^2/(8*Z);                     % power
      P0(3) =  max(ee);                       % field level
    end
  end
else                    % electric walls
  Hx = eta0*interp1(r,H,x,'cubic');
  U  = sum(w.*abs(Hx));
  if sum(H) <= 0
    if gtype == 1
      P0    = -max(hh);                       % field level
    else
      P0(1) = -U;                             % voltage
      P0(2) = -U^2/(8*Z);                     % power
      P0(3) = -max(hh);                       % field level
    end
  else
    if gtype == 1
      P0    =  max(hh);                       % field level
    else
      P0(1) =  U;                             % voltage
      P0(2) =  U^2/(8*Z);                     % power 
      P0(3) =  max(hh);                       % field level
    end
  end
end
% ------------------------------------------------------------  
