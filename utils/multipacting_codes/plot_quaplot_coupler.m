% Function program plot_quaplot_coupler(s)
% --------------------------------------------------------------------
% Plots the counter function, the impact energy, the total counter functions
% and the enhanced counter function for a coupler geometry.
% INPUT  s : 3 - voltage
%            2 - field level
%            1 - power
% --------------------------------------------------------------------
% CALLS TO : load_output_coupler.m, error_message.m
% 09/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% 25/09/00 :                - plot for the total counter added
% ---------------------------------------------------------------------

function plot_quaplot_coupler(s,ss)

if nargin < 1
  s = 1;
end
if nargin < 2
  ss = 1;
end

load_output_coupler

if ok > 0

  if n == 0
    error_message('Unable to plot the counters. No initial points.');
    return;
  end

  if ok1*ok2 == 0
    error_message('Counter functions or impact energy missing.');
  else
    if ss > 0
      cl=error_message(['Plotting the quaplot (counter functions and '...
		       'impact energy).']);
    end

    Window = figure(5);
    clf;
    set(Window,'name','MULTIPAC - Output Window I');

    figure(5)
    if s == 3
      subplot(4,1,1)
      plot(U/1e3,C/n),grid
      ylabel(['c_{' num2str(N) '} / c_0 ' ])
      xlabel('Voltage  [kV]')
      title(['MultiPac 2.1             Counter Function             ' date ])
      ax = axis;
      axis([min(U)/1e3,max(U)/1e3,0,max([0.1,ax(4)])]);

      subplot(4,1,2)
      plot(U/1e3,Efq),grid
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

      subplot(4,1,3)
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

      subplot(4,1,4)
      semilogy(U/1e3,(At+1)/n),grid  
      xlabel('Voltage   [kV]')
      hold on
      plot([min(U)/1e3,max(U)/1e3],[1,1],'-r')
      hold off
      ax = axis;
      axis([min(U)/1e3,max(U)/1e3,min((At+1)/n),ax(4)]);
      ylabel(['t_{' num2str(N) '} / c_0 ' ])
      xlabel('Voltage  [kV]')
      title('Total Counter Function')
    
    
    elseif s == 2
  
      subplot(4,1,1)
      plot(Efl/1e3,C/n),grid
      ylabel(['c_{' num2str(N) '} / c_0 ' ])
      xlabel('Peak electric field  [kV/m]')
      title(['MultiPac 2.1             Counter function             ' date ])
      ax = axis;
      axis([min(Efl)/1e3,max(Efl)/1e3,0,max([0.1,ax(4)])]);

      subplot(4,1,2)
      plot(Efl/1e3,Efq),grid
      hold on      
%      plot([min(Efl)/1e3,max(Efl)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
      e0 = interp1(secy1(1:e1,2),secy1(1:e1,1),1);
      plot([min(Efl)/1e3,max(Efl)/1e3],[e0,e0],'-r')
      plot([min(Efl)/1e3,max(Efl)/1e3],[secy1(e2,1),secy1(e2,1)],'-r')
      plot([min(Efl)/1e3,max(Efl)/1e3],[secy1(e3,1),secy1(e3,1)],'--r')
      hold off
      ylabel(['Ef_{' num2str(N) '} ' ])
      xlabel('Peak electric field  [kV/m]')
      title('Final Impact Energy in eV')
      ax = axis;
      axis([min(Efl)/1e3,max(Efl)/1e3,0,ax(4)]);

      subplot(4,1,3)
      semilogy(Efl/1e3,(A+1)/n),grid  
      xlabel('Voltage   [kV]')
      hold on
      plot([min(Efl)/1e3,max(Efl)/1e3],[1,1],'-r')
      hold off
      ax = axis;
      axis([min(Efl)/1e3,max(Efl)/1e3,min((A+1)/n),ax(4)]);
      ylabel(['e_{' num2str(N) '} / c_0 ' ])
      xlabel('Peak electric field   [kV/m]')
      title('Enhanced Counter Function')

      subplot(4,1,4)
      semilogy(Efl/1e3,(At+1)/n),grid  
      xlabel('Voltage   [kV]')
      hold on
      plot([min(Efl)/1e3,max(Efl)/1e3],[1,1],'-r')
      hold off
      ax = axis;
      axis([min(Efl)/1e3,max(Efl)/1e3,min((At+1)/n),ax(4)]);
      ylabel(['t_{' num2str(N) '} / c_0 ' ])
      xlabel('Peak electric field   [kV/m]')
      title('Total Counter Function')
      
    elseif s == 1

      subplot(4,1,1)
      plot(Pow/1e3,C/n),grid
      ylabel(['c_{' num2str(N) '} / c_0 ' ])
      xlabel('Power  [kW]')
      title(['MultiPac 2.1             Counter function             ' date ])
      ax = axis;
      axis([min(Pow)/1e3,max(Pow)/1e3,0,max([0.1,ax(4)])]);

      subplot(4,1,2)
      plot(Pow/1e3,Efq),grid
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

      subplot(4,1,3)
      semilogy(Pow/1e3,(A+1)/n),grid  
      xlabel('Voltage   [kV]')
      hold on
      plot([min(Pow)/1e3,max(Pow)/1e3],[1,1],'-r')
      hold off
      ax = axis;
      axis([min(Pow)/1e3,max(Pow)/1e3,min((A+1)/n),ax(4)]);
      ylabel(['e_{' num2str(N) '} / c_0 ' ])
      xlabel('Power  [kW]')
      title('Enhanced Counter Function')
    
      subplot(4,1,4)
      semilogy(Pow/1e3,(At+1)/n),grid  
      xlabel('Voltage   [kV]')
      hold on
      plot([min(Pow)/1e3,max(Pow)/1e3],[1,1],'-r')
      hold off
      ax = axis;
      axis([min(Pow)/1e3,max(Pow)/1e3,min((At+1)/n),ax(4)]);
      ylabel(['t_{' num2str(N) '} / c_0 ' ])
      xlabel('Power  [kW]')
      title('Total Counter Function')
    
    end  
  end
end
% --------------------------------------------------------------------
