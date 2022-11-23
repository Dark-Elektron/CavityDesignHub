% Function program run_mpanalysis(s)
% ----------------------------------------------------------------------------
% Carries out the multipacing analysis. Calculates the counter functions, 
% distance map at the maximum of the enhanced counter function and electron 
% trajectory at the minimum of the distance function.
%
% ----------------------------------------------------------------------------
% CALLS TO : clear_window.m, test_geometry.m, check_fieldparam.m,
%            check_inputs.m, check_inputs_win.m, plot_insec.m, error_message.m
%            plot_inisecang.m, mp_cavity_coupler.m, mp_window.m, 
%            check_kama.m
% 09/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% 28/09/00 :                - if initial angle distribution file exists, the
%                             distribution is plotted
% 17/10/00 :                - coefficient files removed
% 13/12/00 :                - check_kama.m added
% ----------------------------------------------------------------------------

function run_mpanalysis(s)

if nargin < 1
  s = 1;
end
if s == 1
  clear_window;
end

ok = test_geometry(0);
if ok == 0
  return;
end

ok0 = check_fieldparam;
if ok0 > 0
  load fieldparam
  gtype = fieldparam(1);
  if gtype <= 2
    ok1 = check_inputs;
  else
    ok1 = check_inputs_win;
  end
else
  ok1 = 0;
end

ok2 = check_kama;

if ok1*ok2 > 0

  figure(1);
  error_message('------------ Multipacting Analysis -----------.');
  error_message('                                               ');
  
  % plot the initial points, the geometry and the secondary yield curve
  if exist('initangle')
    plot_inisecang(1);
  else
    plot_inisec(1);
  end

  % start the MP analysis
  if gtype <= 2
    mp_cavity_coupler(s);
  elseif gtype == 3
    mp_window(s);
  end

  if s > 0
    error_message('                                              ');
    error_message('------ Multipacting analysis completed ------ ');
  end    
end  
% ---------------------------------------------------------------------



