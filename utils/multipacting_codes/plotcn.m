% Function program plotcn
% ---------------------------------------------------------------------------
% Plots the counter function.
%
% ---------------------------------------------------------------------------
% CALLS TO : load_output_data.m, load_output_coupler.m, load_output_window.m
% 06/03/00 : Pasi Yla-Oijala Rolf Nevanlinna Institute
% ---------------------------------------------------------------------------

function plotcn

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
    C = Cl; n = nl;
  elseif side == 'r'
    C = Cr; n = nr;
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
    plot(U/1e3,C),grid
    xlabel('Peak voltage   [kV]')
    ax = axis;
    axis([min(U)/1e3,max(U)/1e3,0,ax(4)]);
  elseif dval == 3
    plot(Efl/1e3,C),grid
    xlabel('Peak Electric Field   [kV/m]')
    ax = axis;
    axis([min(Efl)/1e3,max(Efl)/1e3,0,ax(4)]);
  elseif dval == 1
    plot(Pow/1e3,C),grid
    xlabel('Forward Power   [kW]')
    ax = axis;
    axis([min(Pow)/1e3,max(Pow)/1e3,0,ax(4)]);
  end
  ylabel(['c_{' num2str(N) '} ' ])
  title(['MultiPac 2.0              Counter Function                ' date ])
elseif tval == 2
  if dval == 2
    plot(U/1e3,C/n),grid
    xlabel('Peak voltage   [kV]')
    ax = axis;
    axis([min(U)/1e3,max(U)/1e3,0,ax(4)]);
  elseif dval == 3
    plot(Efl/1e3,C/n),grid
    xlabel('Peak Electric Field   [kV/m]')
    ax = axis;
    axis([min(Efl)/1e3,max(Efl)/1e3,0,ax(4)]);
  elseif dval == 1
    plot(Pow/1e3,C/n),grid
    xlabel('Forward Power   [kW]')
    ax = axis;
    axis([min(Pow)/1e3,max(Pow)/1e3,0,ax(4)]);
  end
  ylabel(['c_{' num2str(N) '} / c_0 ' ])
  title(['MultiPac 2.0          Relative Counter Function             ' date ])
elseif tval == 3
  if dval == 2
    plot(U/1e3,C/n*100),grid
    xlabel('Peak voltage   [kV]')
    ax = axis;
    axis([min(U)/1e3,max(U)/1e3,0,ax(4)]);
  elseif dval == 3
    plot(Efl/1e3,C/n*100),grid
    xlabel('Peak Electric Field   [kV/m]')
    ax = axis;
    axis([min(Efl)/1e3,max(Efl)/1e3,0,ax(4)]);
  elseif dval == 1
    plot(Pow/1e3,C/n*100),grid
    xlabel('Forward Power   [kW]')
    ax = axis;
    axis([min(Pow)/1e3,max(Pow)/1e3,0,ax(4)]);
  end
  ylabel(['c_{' num2str(N) '} / c_0  * 100 %' ])
  title(['MultiPac 2.0        Relative Counter Function in %          ' date ])
end
% --------------------------------------------------------------------



