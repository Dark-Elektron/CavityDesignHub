% Function program plot_initial_angle
% ------------------------------------------------------------------------
% Plots the initial angle distribution (Cumulative distribution fuction).
%
% ------------------------------------------------------------------------
% CALLS TO : error_message.m, clear_window.m
% 25/09/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% ------------------------------------------------------------------------

function plot_initial_angle

cl = 0;
if exist('initangle')
  cl = error_message('Plotting the emission angle distribution.');
  load initangle
  angle1 = initangle(:,1);
  angle2 = initangle(:,2);
  
  Window = figure(3);
  clf;
  set(Window,'name','MULTIPAC - Input Window II');

  subplot(1,1,1)
  plot(angle1*180/pi,angle2,'-r')
  grid
  xlabel('Emission angle [degrees]')
  ylabel('Cumulative distribution fuction')
  title(['MultiPac 2.0                   Emission angle           ' date ])   
else  
  cl = error_message('Emission angle is constant.');
end  

if cl == 1
  clear_window;
end
% ---------------------------------------------------------------------
