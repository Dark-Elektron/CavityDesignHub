% Function program yi0 = plot_distance(trajec,side,s)
% --------------------------------------------------------------------------
% Plots the distance function.
% INPUT  trajec : 1 - determine initial point for the trajectory calculation
%        side   : 'l'eft or 'r'ight, for a window, 'c' otherwise
%
% --------------------------------------------------------------------------
% CALLS TO : clear_window.m, check_inputs.m, check_outputs.m, error_message.m,
%            check_inputs_win.m, check_outputs_win.m, mapplot1.m, 
%            check_geotype.m
% 14/03/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% 13/12/00 :                 - check_geotype.m added
% ---------------------------------------------------------------------------

function yi0 = plot_distance(trajec,side,s)

if nargin < 1
  trajec = 0;
end  
if nargin < 2
  side = 'c';
end
if nargin < 3
  s = 1;
end  

ok4 = check_geotype(side);
if ok4 == 0
  return;
end

ok1 = 1; ok2 = 1; ok3 = 1;
cl = 0;
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

if ok1 ~= 0 
  ok3 = check_fieldparam;
  if ok3 > 0
    load fieldparam
    gtype = fieldparam(1);
  end
end

if (ok1*ok2*ok3) > 0
  if s > 0
    if gtype < 3
      side = 'c';
      cl = error_message('Plotting the distance map.');
    elseif gtype == 3
      if side == 'l' 
        cl = error_message('Plotting the distance map. Warm side.');
      elseif side == 'r'
        cl = error_message('Plotting the distance map. Cold side.');
      else
	cl = error_message(['To plot the distance function, choose warm '...
			    'or cold side.']);
	return;
      end
    end
  end

  Window = figure(6);
  clf;
  set(Window,'name','MULTIPAC - Output Window II');

  load param
  if side == 'l'
    load Ddistancel
    load counter_initialsl
    initials = initialsl;
    load geodatal.n
    ng = length(geodatal(:,1));
    bo = geodatal(4:ng,1:2)';
  elseif side == 'r'
    load Ddistancer
    load counter_initialsr
    initials = initialsr;
    load geodatar.n
    ng = length(geodatar(:,1));
    bo = geodatar(4:ng,1:2)';
  else
    load Ddistance
    load counter_initials
    load geodata.n
    ng = length(geodata(:,1));
    bo = geodata(4:ng,1:2)';
  end

  if min(min(abs(D)==2))
    cl = error_message('Distance map is zero.');
  else
    yi0 = mapplot1(D,1,param,initials,bo,trajec);
  
    N = param(5);
    if side == 'l'
      ti = [' Distance map   d_{' num2str(N) '}   Warm Side'];
      title(['MultiPac 2.0       ' ti '           ' date])    
    elseif side == 'r'
      ti = [' Distance map   d_{' num2str(N) '}   Cold Side'];
      title(['MultiPac 2.0       ' ti '           ' date])
    end
  end
end  

if cl == 1
  clear_window;
end
% --------------------------------------------------------------------
