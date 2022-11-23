% Function program mp_cavity_coupler(s)
% --------------------------------------------------------------------------
% Carries out the multipacing analysis in a cavity and a coupler geometry.
%
% --------------------------------------------------------------------------
% CALLS TO : plot_mixed_fields.m, calculate_counters.m, plot_triplot.m, 
%            plot_triplot_coupler.m, calculate_distance.m, plot_distance.m, 
%            calculate_trajectory.m, plot_trajectory.m, error_message.m
% 09/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% --------------------------------------------------------------------------

function mp_cavity_coupler(s)

if nargin == 0
  s = 1;
end

load fieldparam;
gtype = fieldparam(1);

% plot the electromagnetic fields
if gtype == 2
  plot_mixed_fields(0,'c');
else
  plot_FEM_fields(0,'c');
end  

% calculate the counter functions
calculate_counters(0);

load Ccounter
load Acounter
load counter_flevels

figure(5)
if gtype == 1
  plot_triplot(2);
elseif gtype == 2
  plot_triplot_coupler(1);
  R   = fieldparam(5)+i*fieldparam(6);
  Z   = fieldparam(8);
  Pow = volt2pow(flevel,R,Z);
end  

cl = 0;
if max(C) == 0
  cl=error_message(['Counter function is dentically zero and multipacting ' ...
		    'analysis is completed.']);
else  
  [val,ind] = max(A);
  subplot(3,1,1)
  ax = axis;
  hold on
  if gtype == 1
    plot([flevel(ind)/1e3,flevel(ind)/1e3],[ax(3),ax(4)],'-r')
  elseif gtype == 2
    plot([Pow(ind)/1e3,Pow(ind)/1e3],[ax(3),ax(4)],'-r')
  end  
  hold off

  subplot(3,1,2)
  ax = axis;
  hold on
  if gtype == 1
    plot([flevel(ind)/1e3,flevel(ind)/1e3],[ax(3),ax(4)],'-r')
  elseif gtype == 2
    plot([Pow(ind)/1e3,Pow(ind)/1e3],[ax(3),ax(4)],'-r')
  end  
  hold off

  subplot(3,1,3)
  ax = axis;
  hold on
  if gtype == 1
    plot([flevel(ind)/1e3,flevel(ind)/1e3],[ax(3),ax(4)],'-r')
  elseif gtype == 2
    plot([Pow(ind)/1e3,Pow(ind)/1e3],[ax(3),ax(4)],'-r')
  end  
  hold off

  % calculate the distance map
  calculate_distance(0);

  % plot the distance map
  plot_distance(0,'c',0);
  load Ddistance
  [val,yi0] = min(abs(D));

  if s > 0
    cl = error_message('                                  ');
  end  

  % calculate an electron trajectory
  calculate_trajectory(0,'c');

  % plot the trajectory
  plot_trajectory('c',0);
end  
% ---------------------------------------------------------------------
