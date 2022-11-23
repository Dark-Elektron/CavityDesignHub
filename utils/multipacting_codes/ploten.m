% Function program ploten
% --------------------------------------------------------------------------
% Plots the enhanced counter function.
%
% --------------------------------------------------------------------------
% CALLS TO : load_output_data.m, load_output_coupler.m, load_output_window.m
% 06/03/00 : Pasi Yla-Oijala Rolf Nevanlinna Institute
% --------------------------------------------------------------------------

function ploten

load side

load fieldparam
gtype = fieldparam(1);
if gtype == 1
  load_output_data
elseif gtype == 2
  load_output_coupler
elseif gtype == 3
  load_output_window
  if side == 'l'
    A = Al; n = nl;
  elseif side == 'r'
    A = Ar; n = nr;
  end  
end

htyp = findobj('Tag','Plot type');
if gtype == 1
  dval = 3;
else  
  hdim = findobj('Tag','Dimension');
  dval = get(hdim,'value');
end  
tval = get(htyp,'value');


if tval == 1
  if dval == 2
    plot(U/1e3,A),grid
    xlabel('Peak voltage   [kV]')
    ax = axis;
    axis([min(U)/1e3,max(U)/1e3,ax(3),ax(4)]);
  elseif dval == 3
    plot(Efl/1e3,A),grid
    xlabel('Peak Electric Field   [kV/m]')
    ax = axis;
    axis([min(Efl)/1e3,max(Efl)/1e3,ax(3),ax(4)]);
  elseif dval == 1
    plot(Pow/1e3,A),grid
    xlabel('Forward Power   [kW]')
    ax = axis;
    axis([min(Pow)/1e3,max(Pow)/1e3,ax(3),ax(4)]);
  end
  ylabel(['e_{' num2str(N) '} ' ])
  title(['MultiPac 2.0            Enhanced Counter Function           ' date ])
elseif tval == 2
  if dval == 2
    plot(U/1e3,A/n),grid
    xlabel('Peak voltage   [kV]')
    ax = axis;
    axis([min(U)/1e3,max(U)/1e3,ax(3),ax(4)]);
  elseif dval == 3
    plot(Efl/1e3,A/n),grid
    xlabel('Peak Electric Field   [kV/m]')
    ax = axis;
    axis([min(Efl)/1e3,max(Efl)/1e3,ax(3),ax(4)]);
  elseif dval == 1
    plot(Pow/1e3,A/n),grid
    xlabel('Forward Power   [kW]')
    ax = axis;
    axis([min(Pow)/1e3,max(Pow)/1e3,ax(3),ax(4)]);
  end
  ylabel(['e_{' num2str(N) '} / c_0 ' ])
  title(['MultiPac 2.0       Relative Enhanced Counter Function           ' date ])
elseif tval == 3  
  if dval == 2
    semilogy(U/1e3,A+1),grid  
    xlabel('Peak voltage   [kV]')
    ax = axis;
    axis([min(U)/1e3,max(U)/1e3,1,ax(4)]);
  elseif dval == 3
    semilogy(Efl/1e3,A+1),grid  
    xlabel('Peak Electric Field   [kV/m]')
    ax = axis;
    axis([min(Efl)/1e3,max(Efl)/1e3,1,ax(4)]);
  elseif dval == 1
    semilogy(Pow/1e3,A+1),grid  
    xlabel('Forward Power   [kW]')
    ax = axis;
    axis([min(Pow)/1e3,max(Pow)/1e3,1,ax(4)]);
  end
  ylabel(['e_{' num2str(N) '} ' ])
  title(['MultiPac 2.0       Enhanced Counter Function in log_{10}-scale      ' date ])
elseif tval == 4
  if dval == 2
    semilogy(U/1e3,(A+1)/n),grid  
    xlabel('Peak voltage   [kV]')
    hold on
    plot([min(U)/1e3,max(U)/1e3],[1,1],'-r')
    hold off
    ax = axis;
    axis([min(U)/1e3,max(U)/1e3,min((A+1)/n),ax(4)]);
  elseif dval == 3
    semilogy(Efl/1e3,(A+1)/n),grid  
    xlabel('Peak Electric Field   [kV/m]')
    hold on
    plot([min(Efl)/1e3,max(Efl)/1e3],[1,1],'-r')
    hold off
    ax = axis;
    axis([min(Efl)/1e3,max(Efl)/1e3,min((A+1)/n),ax(4)]);
  elseif dval == 1
    semilogy(Pow/1e3,(A+1)/n),grid  
    xlabel('Forward Power   [kW]')
    hold on
    plot([min(Pow)/1e3,max(Pow)/1e3],[1,1],'-r')
    hold off
    ax = axis;
    axis([min(Pow)/1e3,max(Pow)/1e3,min((A+1)/n),ax(4)]);
  end
  ylabel(['e_{' num2str(N) '} / c_0 ' ])
  title(['MultiPac 2.0     Relative Enhanced Counter Function in log_{10}-scale      ' date ])
end
% --------------------------------------------------------------------
