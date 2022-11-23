% Function program plot_energy(side)
% ---------------------------------------------------------------------------
% Plots the final impact energy.
% INPUT  : 'l'eft or 'r'ight, for a window
% ---------------------------------------------------------------------------
% CALLS TO : clear_window.m, check_inputs.m, check_outputs.m, error_message.m,
%            check_inputs_win.m, check_outputs_win.m, load_output_data.m, 
%            load_output_coupler.m, load_output_window.m, plotef.m,
%            check_geotype.m
% 10/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% 13/12/00 :                - check_geotype.m added
% ---------------------------------------------------------------------------

function plot_energy(side)

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
  load fieldparam
  gtype = fieldparam(1);
  
  if gtype == 1
    load_output_data
    if n == 0
      cl=error_message('Unable to plot the impact energy. No initial points.');
      ok = 0;
    end
    
    if ok > 0
      cl = error_message('Plotting the final impact energy.');

      Window = figure(5);
      clf;
      set(Window,'name','MULTIPAC - Output Window I');

      subplot(1,1,1)
      plot(flevel/1e3,Efq),grid
      hold on
%      plot([min(flevel)/1e3,max(flevel)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
      e0 = interp1(secy1(1:e1,2),secy1(1:e1,1),1);
      plot([min(flevel)/1e3,max(flevel)/1e3],[e0,e0],'-r')
      plot([min(flevel)/1e3,max(flevel)/1e3],[secy1(e2,1),secy1(e2,1)],'-r')
      plot([min(flevel)/1e3,max(flevel)/1e3],[secy1(e3,1),secy1(e3,1)],'--r')

      if exist('secy2')
	load secy2
        e21 = min(find(secy2(:,2)>=1));      % lower thershold
        e20 =  interp1(secy2(1:e21,2),secy2(1:e21,1),1);
        e22 = max(find(secy2(:,2)>=1));      % upper thershold
        [val,e23] = max(secy2(:,2));         % maximum secondary yield
	plot([min(flevel)/1e3,max(flevel)/1e3],[e20,e20],'-b')
       plot([min(flevel)/1e3,max(flevel)/1e3],[secy2(e22,1),secy2(e22,1)],'-b')
      plot([min(flevel)/1e3,max(flevel)/1e3],[secy2(e23,1),secy2(e23,1)],'--b')
      end
      if exist('secy3')
	load secy3
        e31 = min(find(secy3(:,2)>=1));      % lower thershold
        e30 =  interp1(secy3(1:e31,2),secy3(1:e31,1),1);
        e32 = max(find(secy3(:,2)>=1));      % upper thershold
        [val,e33] = max(secy3(:,2));         % maximum secondary yield
	plot([min(flevel)/1e3,max(flevel)/1e3],[e30,e30],'-g')
       plot([min(flevel)/1e3,max(flevel)/1e3],[secy3(e32,1),secy3(e32,1)],'-g')
      plot([min(flevel)/1e3,max(flevel)/1e3],[secy3(e33,1),secy3(e33,1)],'--g')
      end
      if exist('secy4')
	load secy4
        e41 = min(find(secy4(:,2)>=1));      % lower thershold
        e40 =  interp1(secy4(1:e41,2),secy1(1:e41,1),1);
        e42 = max(find(secy4(:,2)>=1));      % upper thershold
        [val,e43] = max(secy4(:,2));         % maximum secondary yield
	plot([min(flevel)/1e3,max(flevel)/1e3],[e40,e40],'-k')
       plot([min(flevel)/1e3,max(flevel)/1e3],[secy4(e42,1),secy4(e42,1)],'-k')
      plot([min(flevel)/1e3,max(flevel)/1e3],[secy4(e43,1),secy4(e43,1)],'--k')
      end
      hold off
      ylabel(['Ef_{' num2str(N) '} ' ])
      xlabel('Peak Electric Field  [kV/m]')
      title(['MultiPac 2.1           Final Impact Energy in ev        ' date ])
 
      htyp = uicontrol('Style','Popup','String','lin|sqrt|loq','Units',...
		       'Normalized','Position',[0.78 0.85 0.12 0.05],...
		       'Callback','plotef','Tag','Plot type');

%     hdim = uicontrol('Style','Popup','String','field|voltage|power',...
%	       'Units','Normalized','Position',[0.78 0.78 0.12 0.05],...
%	       'Callback','plotef','Tag','Dimension');
    end
  
  elseif gtype == 2
    load_output_coupler
    if n == 0
      cl=error_message('Unable to plot the impact energy. No initial points.');
      ok = 0;
    end
    
    if ok > 0
      cl = error_message('Plotting the final impact energy.');

      Window = figure(5);
      clf;
      set(Window,'name','MULTIPAC - Output Window I');

      subplot(1,1,1)
      plot(Pow/1e3,Efq),grid
      hold on
%      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
      e0 = interp1(secy1(1:e1,2),secy1(1:e1,1),1);
      plot([min(Pow)/1e3,max(Pow)/1e3],[e0,e0],'-r')
      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e2,1),secy1(e2,1)],'-r')
      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e3,1),secy1(e3,1)],'--r')
      
      if exist('secy2')
	load secy2
        e21 = min(find(secy2(:,2)>=1));      % lower thershold
        e20 =  interp1(secy2(1:e21,2),secy2(1:e21,1),1);
        e22 = max(find(secy2(:,2)>=1));      % upper thershold
        [val,e23] = max(secy2(:,2));         % maximum secondary yield
	plot([min(Pow)/1e3,max(Pow)/1e3],[e20,e20],'-k')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy2(e22,1),secy2(e22,1)],'-k')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy2(e23,1),secy2(e23,1)],'--k')
      end
      if exist('secy3')
	load secy3
        e31 = min(find(secy3(:,2)>=1));      % lower thershold
        e30 =  interp1(secy3(1:e31,2),secy3(1:e31,1),1);
        e32 = max(find(secy3(:,2)>=1));      % upper thershold
        [val,e33] = max(secy3(:,2));         % maximum secondary yield
	plot([min(Pow)/1e3,max(Pow)/1e3],[e30,e30],'-g')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy3(e32,1),secy3(e32,1)],'-g')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy3(e33,1),secy3(e33,1)],'--g')
      end
      if exist('secy4')
	load secy4
        e41 = min(find(secy4(:,2)>=1));      % lower thershold
        e40 =  interp1(secy4(1:e41,2),secy1(1:e41,1),1);
        e42 = max(find(secy4(:,2)>=1));      % upper thershold
        [val,e43] = max(secy4(:,2));         % maximum secondary yield
	plot([min(Pow)/1e3,max(Pow)/1e3],[e40,e40],'-m')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy4(e42,1),secy4(e42,1)],'-m')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy4(e43,1),secy4(e43,1)],'--m')
      end
      
      hold off
      ylabel(['Ef_{' num2str(N) '} ' ])
      xlabel('RF power  [kW]')
      title(['MultiPac 2.1          Final Impact Energy in eV         ' date ])
 
      htyp = uicontrol('Style','Popup','String','lin|sqrt|loq','Units',...
		       'Normalized','Position',[0.78 0.85 0.12 0.05],...
		       'Callback','plotef','Tag','Plot type');

      hdim = uicontrol('Style','Popup','String','power|voltage',...
		      'Units','Normalized','Position',[0.78 0.78 0.12 0.05],...
		      'Callback','plotef','Tag','Dimension');
    end
    
  elseif gtype == 3
    load_output_window
    
    if ok > 0
      if side == 'l'
	Efq = Eql; n = nl;
      elseif side == 'r'
	Efq = Eqr; n = nr;
      end
      if n == 0
	cl = error_message(['Unable to plot the impact energy. No ' ...
			    'initial points.']);
	ok = 0;
	return;
      end
      
      if side == 'l'
        cl = error_message('Plotting the final impact energy. Warm side.');
      elseif side == 'r'
        cl = error_message('Plotting the final impact energy. Cold side.'); 
      else
        cl = error_message(['To plot the impact energy in window, ' ...
			    'choose warm or cold side.']);
	return;
      end

      Window = figure(5);
      clf;
      set(Window,'name','MULTIPAC - Output Window I');

      subplot(1,1,1)
      plot(Pow/1e3,Efq),grid
      hold on
%      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
      e0 = interp1(secy1(1:e1,2),secy1(1:e1,1),1);
      plot([min(Pow)/1e3,max(Pow)/1e3],[e0,e0],'-r')
      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e2,1),secy1(e2,1)],'-r')
      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e3,1),secy1(e3,1)],'--r')

      if exist('secy2')
	load secy2
        e21 = min(find(secy2(:,2)>=1));      % lower thershold
        e20 =  interp1(secy2(1:e21,2),secy2(1:e21,1),1);
        e22 = max(find(secy2(:,2)>=1));      % upper thershold
        [val,e23] = max(secy2(:,2));         % maximum secondary yield
	plot([min(Pow)/1e3,max(Pow)/1e3],[e20,e20],'-k')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy2(e22,1),secy2(e22,1)],'-k')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy2(e23,1),secy2(e23,1)],'--k')
      end
      if exist('secy3')
	load secy3
        e31 = min(find(secy3(:,2)>=1));      % lower thershold
        e30 =  interp1(secy3(1:e31,2),secy3(1:e31,1),1);
        e32 = max(find(secy3(:,2)>=1));      % upper thershold
        [val,e33] = max(secy3(:,2));         % maximum secondary yield
	plot([min(Pow)/1e3,max(Pow)/1e3],[e30,e30],'-g')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy3(e32,1),secy3(e32,1)],'-g')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy3(e33,1),secy3(e33,1)],'--g')
      end
      if exist('secy4')
	load secy4
        e41 = min(find(secy4(:,2)>=1));      % lower thershold
        e40 =  interp1(secy4(1:e41,2),secy1(1:e41,1),1);
        e42 = max(find(secy4(:,2)>=1));      % upper thershold
        [val,e43] = max(secy4(:,2));         % maximum secondary yield
	plot([min(Pow)/1e3,max(Pow)/1e3],[e40,e40],'-m')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy4(e42,1),secy4(e42,1)],'-m')
	plot([min(Pow)/1e3,max(Pow)/1e3],[secy4(e43,1),secy4(e43,1)],'--m')
      end
            
      hold off
      ax = axis;
      axis([min(Pow)/1e3,max(Pow)/1e3,0,ax(4)]);
      ylabel(['Ef_{' num2str(N) '} ' ])
      xlabel('RF power  [kW]')
      title(['MultiPac 2.1           Final Impact Energy in eV        ' date ])
 
      htyp = uicontrol('Style','Popup','String','lin|sqrt|loq','Units',...
		       'Normalized','Position',[0.78 0.85 0.12 0.05],...
		       'Callback','plotef','Tag','Plot type');

      hdim = uicontrol('Style','Popup','String','power|voltage',...
	       'Units','Normalized','Position',[0.78 0.78 0.12 0.05],...
		       'Callback','plotef','Tag','Dimension');
    end
  end
end  

if cl == 1
  clear_window;
end
% --------------------------------------------------------------------
