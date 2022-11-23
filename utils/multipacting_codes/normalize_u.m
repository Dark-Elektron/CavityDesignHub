% Function program normalize_u(wall)
% ---------------------------------------------------------------------
% Saves and normalizes the calculated fields. 
% INPUT  wall : 0 - cavity
%               1 - coupler of window with E-walls
%               2 - coupler of window with H-walls
% ---------------------------------------------------------------------
% CALLS TO : peak_cavity_field.m, peak_coupler_field.m
% 10/10/00 : Jani Lukkarinen, Rolf Nevanlinna Institute
% ---------------------------------------------------------------------

function normalize_u(wall)

load job
if wall == 0
  E0 = peak_cavity_field
  % normalize u
  load kama0.mat
  u = u ./ E0(1);
  save kama_n.mat index -v4
  save kama_n.mat u -append
  save kama_n.mat k -append
%   save kama_n.mat u -v4 -append
%   save kama_n.mat k -v4 -append
elseif wall == 1
  E0 = peak_coupler_field(1)
  load kama1.mat
  u = u ./ E0(1);
  save kama_1.mat index -v4
  save kama_1.mat u -append
  save kama_1.mat k -append
%   save kama_1.mat u -v4 -append
%   save kama_1.mat k -v4 -append
else 
  E0 = peak_coupler_field(2)
  load kama2.mat
  u = u ./ E0(1);
  save kama_2.mat index -v4
  save kama_2.mat u -append
  save kama_2.mat k -append 
%   save kama_2.mat u -v4 -append
%   save kama_2.mat k -v4 -append  
end
% ------------------------------------------------------------------------
