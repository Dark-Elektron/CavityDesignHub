% Function program run_field_solver.m
% ----------------------------------------------------------------------------
% Runs the field solver.
%
% ----------------------------------------------------------------------------
% CALLS TO : clear_window.m, test_geometry.m, error_message.m, 
%            check_geometry.m, cavity_field.m, coupler_field.m, window_field.m
% 04/04/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% ----------------------------------------------------------------------------

function run_field_solver(s)

if nargin < 1
  s = 1;  
end
if s == 1
  clear_window;
end  

okt = test_geometry(0);
if okt == 0
  return;
end

ok = check_geometry;

if ok > 0
  figure(1);
  error_message('---------- Field Solver ----------');
  error_message('                                  ');

  load fieldparam
  gtype = fieldparam(1);

  if gtype == 1
    cavity_field(s);
  elseif gtype == 2
    coupler_field(s);
  elseif gtype == 3
    window_field(s);
  end  

  if s > 0
%    error_message('                           ');    
    error_message('--------- The end ---------');
  end
  error_message('                          ');
end  
% ---------------------------------------------------------------------

