% Function program plot_quaplot(s)
% -----------------------------------------------------------------------
% Plots the counter function, the impact energy, total counter functions
% and the enhanced counter function.
%
% -----------------------------------------------------------------------
% CALLS TO : load_output_data.m, error_message.m
% 06/03/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% 25/09/00 :                 - plot for total counter added
% ------------------------------------------------------------------------

function plot_quaplot(ss)

if nargin < 1
  ss = 1;
end

load_output_data

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
%    plot([min(Efl)/1e3,max(Efl)/1e3],[secy1(e1,1),secy1(e1,1)],'-r')
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
    title('Enhanced counter function')  
  
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
    title('Total counter function')  
  
  end
end
% --------------------------------------------------------------------
