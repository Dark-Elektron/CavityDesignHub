% Function program plot_triplot_win(s,side,ss)
% --------------------------------------------------------------------
% Plots the counter function, the impact energy and the enhanced
% counter function for a window geometry.
% INPUT  s : 3 - voltage
%            2 - field level
%            1 - power
%        side : 'l'eft or 'r'ight for a window
% --------------------------------------------------------------------
% CALLS TO : load_output_window.m, error_message.m, clear_window.m
% 09/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% ---------------------------------------------------------------------

function plot_triplot_win(s,side,ss)

if nargin < 1
  s = 1;
end
if nargin < 2
  load side
elseif nargin == 2
  save side side
end  
if nargin < 3
  ss = 1;
end

cl = 0;
load_output_window

if ok > 0

  if side == 'l'
    C = Cl; A = Al; Eq = Eql;
    n = nl; 
  elseif side == 'r'
    C = Cr; A = Ar; Eq = Eqr;
    n = nr;
  end
  if n == 0
    error_message('Unable to plot the triplot. No initial points.');
    return;
  end

  if ok1*ok2 == 0
    cl = error_message('Counter functions or impact energy missing.');
  else  

    if ss > 0 & side == 'l'
      cl = error_message(['Plotting the triplot (counter, enhanced ' ...
			  'counter and impact energy). Warm side.']);      
    elseif ss > 0 & side == 'r'
      cl = error_message(['Plotting the triplot (counter, enhanced ' ...
			  'counter and impact energy). Cold side.']);
    end

    Window = figure(5);
    clf;
    set(Window,'name','MULTIPAC - Output Window I');

    figure(5)
    if s == 3
      subplot(3,1,1)
      plot(U/1e3,C/n),grid
      ylabel(['c_{' num2str(N) '} / c_0 ' ])
      xlabel('Voltage  [kV]')
      if side == 'l'
	title(['MultiPac 2.1       Counter Function   Warm Side       ' date ])
      elseif side == 'r'
	title(['MultiPac 2.1       Counter Function   Cold Side       ' date ])
      end
      ax = axis;
      axis([min(U)/1e3,max(U)/1e3,0,max([0.1,ax(4)])]);

      subplot(3,1,2)
      plot(U/1e3,Eq),grid
      hold on
%      plot([min(flevel)/1e3,max(flevel)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
      e0 = interp1(secy1(1:e1,2),secy1(1:e1,1),1);      
      plot([min(flevel)/1e3,max(flevel)/1e3],[e0,e0],'-r')
      plot([min(flevel)/1e3,max(flevel)/1e3],[secy1(e2,1),secy1(e2,1)],'-r')
      plot([min(flevel)/1e3,max(flevel)/1e3],[secy1(e3,1),secy1(e3,1)],'--r')
      hold off
      ylabel(['Ef_{' num2str(N) '} ' ])
      xlabel('Voltage [kV]')
      title('Final Impact Energy in eV')
      ax = axis;
      axis([min(U)/1e3,max(U)/1e3,0,ax(4)]);

      subplot(3,1,3)
      semilogy(U/1e3,(A+1)/n),grid  
      xlabel('Voltage   [kV]')
      hold on
      plot([min(U)/1e3,max(U)/1e3],[1,1],'-r')
      hold off
      ax = axis;
      axis([min(U)/1e3,max(U)/1e3,min((A+1)/n),ax(4)]);
      ylabel(['e_{' num2str(N) '} / c_0 ' ])
      xlabel('Voltage  [kV]')
      title('Enhanced Counter Function')
    
    elseif s == 2
  
      subplot(3,1,1)
      plot(Efl/1e6,C/n),grid
      ylabel(['c_{' num2str(N) '} / c_0 ' ])
      xlabel('Peak electric field  [MV/m]')
      if side == 'l'
	title(['MultiPac 2.1       Counter Function   Warm Side       ' date ])
      elseif side == 'r'
	title(['MultiPac 2.1       Counter Function   Cold Side       ' date ])
      end
      ax = axis;
      axis([min(Efl)/1e6,max(Efl)/1e6,0,max([0.1,ax(4)])]);

      subplot(3,1,2)
      plot(Efl/1e6,Eq),grid
      hold on
%      plot([min(Efl)/1e3,max(Efl)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
      e0 = interp1(secy1(1:e1,2),secy1(1:e1,1),1);
      plot([min(Efl)/1e6,max(Efl)/1e6],[e0,e0],'-r')
      plot([min(Efl)/1e6,max(Efl)/1e6],[secy1(e2,1),secy1(e2,1)],'-r')
      plot([min(Efl)/1e6,max(Efl)/1e6],[secy1(e3,1),secy1(e3,1)],'--r')
      hold off
      ylabel(['Ef_{' num2str(N) '} ' ])
      xlabel('Peak electric field  [MV/m]')
      title('Final Impact Energy in eV')
      ax = axis;
      axis([min(Efl)/1e3,max(Efl)/1e3,0,ax(4)]);

      subplot(3,1,3)
      semilogy(Efl/1e6,(A+1)/n),grid  
      xlabel('Voltage   [MV]')
      hold on
      plot([min(Efl)/1e6,max(Efl)/1e6],[1,1],'-r')
      hold off
      ax = axis;
      axis([min(Efl)/1e6,max(Efl)/1e6,min((A+1)/n),ax(4)]);
      ylabel(['e_{' num2str(N) '} / c_0 ' ])
      xlabel('Peak electric field   [MV/m]')
      title('Enhanced counter function')

    elseif s == 1

      subplot(3,1,1)
      plot(Pow/1e3,C/n),grid
      ylabel(['c_{' num2str(N) '} / c_0 ' ])
      xlabel('Power  [kW]')
      if side == 'l'
	title(['MultiPac 2.1       Counter Function   Warm Side       ' date ])
      elseif side == 'r'
	title(['MultiPac 2.1       Counter Function   Cold Side       ' date ])
      end
      ax = axis;
      axis([min(Pow)/1e3,max(Pow)/1e3,0,max([0.1,ax(4)])]);

      subplot(3,1,2)
      plot(Pow/1e3,Eq),grid
      hold on
%      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
      e0 = interp1(secy1(1:e1,2),secy1(1:e1,1),1);      
      plot([min(Pow)/1e3,max(Pow)/1e3],[e0,e0],'-r')
      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e2,1),secy1(e2,1)],'-r')
      plot([min(Pow)/1e3,max(Pow)/1e3],[secy1(e3,1),secy1(e3,1)],'--r')
      hold off
      ylabel(['Ef_{' num2str(N) '} ' ])
      xlabel('Power [kW]')
      title('Final Impact Energy in eV')
      ax = axis;
      axis([min(Pow)/1e3,max(Pow)/1e3,0,ax(4)]);

      subplot(3,1,3)
      semilogy(Pow/1e3,(A+1)/n),grid  
      xlabel('Voltage   [kV]')
      hold on
      plot([min(Pow)/1e3,max(Pow)/1e3],[1,1],'-r')
      hold off
      ax = axis;
      axis([min(Pow)/1e3,max(Pow)/1e3,min((A+1)/n),ax(4)]);
      ylabel(['e_{' num2str(N) '} / c_0 ' ])
      xlabel('Power  [kW]')
      title('Enhanced counter function')
    end  
  end
end

if cl == 1
  clear_window;
end
% --------------------------------------------------------------------
