% Function program plot_enhanced(side)
% --------------------------------------------------------------------
% Plots the enhanced counter function.
% INPUT  : 'l'eft or 'r'ight, for a window
%
% --------------------------------------------------------------------
% CALLS TO : check_inputs.m, check_outputs.m, check_inputs_win.m, 
%            check_outputs_win.m, load_output_data.m, error_message.m,
%            load_output_coupler.m, load_output_window.m, ploten.m
%            check_geotype.m
% 10/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% 13/12/00 :                 - check_geotype.m added
% ---------------------------------------------------------------------

function plot_enhanced(side)

if nargin < 1
  if exist('side')
    load side
  else
    side = 'c';
  end
end
save side side

ok3 = check_geotype(side);
if ok3 == 0
  return;
end

if side == 'l' | side == 'r'
  ok1 = check_inputs_win;
  ok2 = check_outputs_win;
else
  ok1 = check_inputs;
  ok2 = check_outputs;
end

cl = 0;
if ok1*ok2 > 0

  load fieldparam;
  gtype = fieldparam(1);
  
  if gtype == 1
    load_output_data
    if n == 0
      cl = error_message(['Unable to plot the enhanced counter. No ' ...
			  'initial points.']);
      ok = 0;
    end
    
    if ok > 0
      cl = error_message('Plotting the enhanced counter function.');
    
      Window = figure(5);
      clf;
      set(Window,'name','MULTIPAC - Output Window I');

      subplot(1,1,1)
      plot(flevel/1e3,A),grid
      ylabel(['e_{' num2str(N) '} ' ])
      xlabel('Peak Electric Field  [kV/m]')
      title(['MultiPac 2.1          Enhanced Counter Function         ' date ])
      ax = axis;
      axis([min(flevel)/1e3,max(flevel)/1e3,ax(3),ax(4)]);

      htyp = uicontrol('Style','Popup','String','lin|rlin|log|rlog','Units',...
		       'Normalized','Position',[0.78 0.85 0.12 0.05],...
		       'Callback','ploten','Tag','Plot type');

%      hdim = uicontrol('Style','Popup','String','field|voltage|power',...
%	       'Units','Normalized','Position',[0.78 0.78 0.12 0.05],...
%		       'Callback','ploten','Tag','Dimension');
    end
  elseif gtype == 2
    load_output_coupler
    if n == 0
      cl = error_message(['Unable to plot the enhanced counter. No ' ...
			  'initial points.']);
      ok = 0;
    end
    
    if ok > 0
      cl = error_message('Plotting the enhanced counter function.');
    
      Window = figure(5);
      clf;
      set(Window,'name','MULTIPAC - Output Window I');

      subplot(1,1,1)
      plot(Pow/1e3,A),grid
      ylabel(['e_{' num2str(N) '} ' ])
      xlabel('RF power  [kW]')
      title(['MultiPac 2.1          Enhanced Counter Function         ' date ])
      ax = axis;
      axis([min(Pow)/1e3,max(Pow)/1e3,ax(3),ax(4)]);

      htyp = uicontrol('Style','Popup','String','lin|rlin|log|rlog','Units',...
		       'Normalized','Position',[0.78 0.85 0.12 0.05],...
		       'Callback','ploten','Tag','Plot type');

      hdim = uicontrol('Style','Popup','String','power|voltage',...
	       'Units','Normalized','Position',[0.78 0.78 0.12 0.05],...
		       'Callback','ploten','Tag','Dimension');
    end
  elseif gtype == 3
    load_output_window
    
    if ok > 0
      if side == 'l'
	A = Al; n = nl;
      elseif side == 'r'
	A = Ar; n = nr;
      end
      if n == 0
	cl = error_message(['Unable to plot the enhanced counter. No ' ...
			    'initial points.']);
	ok = 0;
	return;
      end

      if side == 'l'
        cl=error_message('Plotting the enhanced counter function. Warm side.');
      elseif side == 'r'
        cl=error_message('Plotting the enhanced counter function. Cold side.');
      else
        cl=error_message(['To plot the enhanced counter function in ' ...
			  'window, choose warm or cold side.']);
	return;
      end
    
      Window = figure(5);
      clf;
      set(Window,'name','MULTIPAC - Output Window I');

      subplot(1,1,1)
      plot(Pow/1e3,A),grid
      ylabel(['e_{' num2str(N) '} ' ])
      xlabel('RF power  [kW]')
      title(['MultiPac 2.1          Enhanced Counter Function         ' date ])
      ax = axis;
      axis([min(Pow)/1e3,max(Pow)/1e3,ax(3),ax(4)]);

      htyp = uicontrol('Style','Popup','String','lin|rlin|log|rlog','Units',...
		       'Normalized','Position',[0.78 0.85 0.12 0.05],...
		       'Callback','ploten','Tag','Plot type');

      hdim = uicontrol('Style','Popup','String','power|voltage',...
	       'Units','Normalized','Position',[0.78 0.78 0.12 0.05],...
		       'Callback','ploten','Tag','Dimension');
    end
  end
end  

if cl == 1
  clear_window;
end
% --------------------------------------------------------------------
