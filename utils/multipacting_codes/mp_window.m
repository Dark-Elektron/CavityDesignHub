% Function program mp_window(s)
% -------------------------------------------------------------------------
% Carries out the multipacing analysis in a window geometry.
%
% -------------------------------------------------------------------------
% CALLS TO : plot_mixed_fields.m, error_message.m, calculate_counters.m, 
%            plot_triplot_win.m, calculate_distance_win.m, plot_distance.m, 
%            calculate_trajectory.m, plot_trajectory.m
% 14/03/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% -------------------------------------------------------------------------

function mp_window(s)

if nargin < 1
  s = 0;
end

% plot the electromagnetic fields
plot_mixed_fields(0,'l');

% calculate the counter functions (for both sides)
calculate_counters_win(0);

if s > 0
  error_message('Triplot for the warm side.');
end  
plot_triplot_win(1,'l');

load Acounterl
load Ccounterl
if max(C) == 0
  error_message('Counter function on the warm side is identically zero.');
else
  load counter_flevels
  [val,ind] = max(A);
  load fieldparam
  Pow = volt2pow(flevel,fieldparam(5)+i*fieldparam(6),fieldparam(8));

  Window = figure(5);
  set(Window,'name','MULTIPAC - Output Window I');
  figure(5)
  
  subplot(3,1,1)
  ax = axis;
  hold on
  plot([Pow(ind)/1e3,Pow(ind)/1e3],[ax(3),ax(4)],'-r')
  hold off

  subplot(3,1,2)
  ax = axis;
  hold on
  plot([Pow(ind)/1e3,Pow(ind)/1e3],[ax(3),ax(4)],'-r')
  hold off

  subplot(3,1,3)
  ax = axis;
  hold on
  plot([Pow(ind)/1e3,Pow(ind)/1e3],[ax(3),ax(4)],'-r')
  hold off

  % calculate the distance map
  calculate_distance_win(0,'l');

  % plot the distance map
  plot_distance(0,'l',0);
  load Ddistancel
  [val,yi0] = min(abs(D));

  if s > 0
    error_message('                                  ');
  end  

  % calculate an electron trajectory
  calculate_trajectory(0,'l');

  % plot the trajectory
  plot_trajectory('l',0);
end  

if s > 0
  error_message('                                    ');
end  
error_message('To see the cold side, press any key.');
pause

plot_triplot_win(1,'r');
if s > 0
  error_message('Triplot for the cold side.');
end  
  
load Acounterr
load Ccounterr
if max(C) == 0
  error_message('Counter function on the cold side is identically zero.');
else
  load counter_flevels
  [val,ind] = max(A);

  load fieldparam
  Pow = volt2pow(flevel,fieldparam(5)+i*fieldparam(6),fieldparam(8));

  Window = figure(5);
  set(Window,'name','MULTIPAC - Output Window I');
  figure(5)
  
  subplot(3,1,1)
  ax = axis;
  hold on
  plot([Pow(ind)/1e3,Pow(ind)/1e3],[ax(3),ax(4)],'-r')
  hold off

  subplot(3,1,2)
  ax = axis;
  hold on
  plot([Pow(ind)/1e3,Pow(ind)/1e3],[ax(3),ax(4)],'-r')
  hold off

  subplot(3,1,3)
  ax = axis;
  hold on
  plot([Pow(ind)/1e3,Pow(ind)/1e3],[ax(3),ax(4)],'-r')
  hold off

  % calculate the distance map
  calculate_distance_win(0,'r');

  % plot the distance map
  plot_distance(0,'r',0);
  load Ddistancer
  [val,yi0] = min(abs(D));

  if s > 0
    error_message('                                  ');
  end  

  % calculate an electron trajectory
  calculate_trajectory(0,'r');

  % plot the trajectory
  plot_trajectory('r',0);
end  
% ---------------------------------------------------------------------
