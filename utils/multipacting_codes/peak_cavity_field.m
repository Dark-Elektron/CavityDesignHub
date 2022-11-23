% Function program E0 = peak_cavity_field
% ------------------------------------------------------------------
% Computes the peak electric field on the surface of a cavity.
%
% OUTPUT  E0 : peak electric field, [V/M]
% ------------------------------------------------------------------
% CALLS TO : field_bin
% 05/04/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
% ------------------------------------------------------------------

function E0 = peak_cavity_field

% load mesh and field solution
load mesh
load kama0

% compute the field on the boundary
load geodata.n

n   = length(geodata(:,1));
ind = find(geodata(4:n,3)==1); 
ind = ind(2:length(ind));
r   = geodata(4:n,1); %r = r(ind)-5e-4; 
r   = r(ind)-1.75e-4;                    % move points inside
z   = geodata(4:n,2); z = z(ind);

rind = 1:length(r);
zind = 1:length(z);

alue = 0;

save zr z -v4
save zr r -append
save zr alue -append
save zr rind -append
save zr zind -append
% save zr r -v4 -append
% save zr alue -v4 -append
% save zr rind -v4 -append
% save zr zind -v4 -append

% compute the field at generated points
!copy /y kama0.mat kama.mat
!Multipac fields

load Er
load Ez
ee = sqrt(abs(Er).^2+abs(Ez).^2);
E0 = max(ee);
% ------------------------------------------------------------  
