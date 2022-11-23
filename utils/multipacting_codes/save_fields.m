% Function program save_fields(wall)
% ---------------------------------------------------------------------
% Saves and normalizes the calculated fields. The fields are saved to
% files fieldfile1.n or fieldfileE.n and fieldfileH.n
% INPUT  wall : 0 - cavity
%               1 - E-walls
%               2 - H-walls
% ---------------------------------------------------------------------
% CALLS TO : peak_cavity_field.m, peak_coupler_field.m
% 09/05/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% ---------------------------------------------------------------------

function save_fields(wall)

load fieldfile1.txt
file = fieldfile1;

% compute the peak electric field on the boundary or the rf power
load job
if wall == 0 & job == 1
  E0 = peak_cavity_field;
else
  E0 = peak_coupler_field(wall);
end  

% normalize the fields
my0 = 4e-7*pi;
er  = find(file(:,1)==0);
ez  = find(file(:,1)==1);
bp  = find(file(:,1)==2);
file(er,6:7) = file(er,6:7) / E0(1);
file(ez,6:7) = file(ez,6:7) / E0(1);
file(bp,6:7) = my0 * file(bp,6:7) / E0(1);
fieldfile1   = file;

E0

if wall == 0
  save -ascii fieldfile1.n fieldfile1
elseif wall == 1
  save -ascii fieldfileE.n fieldfile1
elseif wall == 2
  save -ascii fieldfileH.n fieldfile1
end
% ------------------------------------------------------------------------
