% Function program plot_absorbed_energy(side,s)
% --------------------------------------------------------------------------
% Plots the absorbed energy.
% INPUT  side   : 'l'eft or 'r'ight, for a window, 'c' otherwise
%
% --------------------------------------------------------------------------
% CALLS TO : clear_window.m, check_inputs.m, check_outputs.m, error_message.m,
%            check_inputs_win.m, check_outputs_win.m, mapplot1.m
% 25/09/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% ---------------------------------------------------------------------------

function plot_absorbed_energy(side,s)

if nargin < 1
  side = 'c';
end
if nargin < 2
  s = 1;
end  

if side == 'l'
  ok2 = check_distance_win('l');
  ok1 = check_inputs_win;
elseif side == 'r'
  ok2 = check_distance_win('r');
  ok1 = check_inputs_win;
else
  ok2 = check_distance;
  ok1 = check_inputs;
end

if (ok1*ok2) > 0
  if s > 0
    error_message('Plotting the absorbed energy map.');
  end

  Window = figure(6);
  clf;
  set(Window,'name','MULTIPAC - Output Window II');

  load param
  if side == 'l'
    load Edistancel
    load counter_initialsl
    initials = initialsl;
    load geodatal.n
    ng = length(geodatal(:,1));
    bo = geodatal(4:ng,1:2)';
  elseif side == 'r'
    load Edistancer
    load counter_initialsr
    initials = initialsr;
    load geodatar.n
    ng = length(geodatar(:,1));
    bo = geodatar(4:ng,1:2)';
  else
    load Edistance
    load counter_initials
    load geodata.n
    ng = length(geodata(:,1));
    bo = geodata(4:ng,1:2)';
  end

  if max(max(abs(Em)==0))
    error_message('Absorbed energy map is zero');
  else  
    yi0 = mapplot1(Em,1,param,initials,bo,0);
  
    N = param(5);
    if side == 'l'
      ti = [' Absorbed Energy   em_{' num2str(N) '}   Warm Side'];
      title(['MultiPac 2.0       ' ti '           ' date])    
    elseif side == 'r'
      ti = [' Absorbed Energy   em_{' num2str(N) '}   Cold Side'];
      title(['MultiPac 2.0       ' ti '           ' date])
    end
  end
end  
% --------------------------------------------------------------------

