% Function program plotef.m
% --------------------------------------------------------------------------
% Plots the final impact energy.
%
% --------------------------------------------------------------------------
% CALLS TO : load_output_data.m, load_output_coupler.m, load_output_window.m
% 06/03/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% --------------------------------------------------------------------------

function plotef(side)

if nargin < 1
  load side
end  

load fieldparam
gtype = fieldparam(1);
if gtype == 1
  load_output_data
elseif gtype == 2
  load_output_coupler
elseif gtype == 3
  load_output_window
  if side == 'l'
    Efq = Eql; n = nl;
  elseif side == 'r'
    Efq = Eqr; n = nr;
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
    plot(U/1e3,Efq),grid
    xlabel('Peak voltage   [kV]')
    ax = axis;
    axis([min(U)/1e3,max(U)/1e3,0,ax(4)]);
  elseif dval == 3
    plot(Efl/1e3,Efq),grid
    xlabel('Peak Electric Field   [kV/m]')
    ax = axis;
    axis([min(Efl)/1e3,max(Efl)/1e3,0,ax(4)]);    
  elseif dval == 1
    plot(Pow/1e3,Efq),grid
    xlabel('Forward Power   [kW]')
    ax = axis;
    axis([min(Pow)/1e3,max(Pow)/1e3,0,ax(4)]);    
  end
  ylabel(['Ef_{' num2str(N) '} ' ])
  title(['MultiPac 2.0            Final Impact Energy in ev          ' date ])  
elseif tval == 2
  if dval == 2
    plot(U/1e3,sqrt(Efq)),grid  
    xlabel('Peak voltage   [kV]')
    ax = axis;
    axis([min(U)/1e3,max(U)/1e3,0,ax(4)]);    
  elseif dval == 3
    plot(Efl/1e3,sqrt(Efq)),grid  
    xlabel('Peak Electric Field   [kV/m]')
    ax = axis;
    axis([min(Efl)/1e3,max(Efl)/1e3,0,ax(4)]);    
  elseif dval == 1
    plot(Pow/1e3,sqrt(Efq)),grid  
    xlabel('Forward Power   [kW]')
    ax = axis;
    axis([min(Pow)/1e3,max(Pow)/1e3,0,ax(4)]);    
  end
  ylabel([ 'sqrt of Ef_{' num2str(N) '} ' ])
  title(['MultiPac 2.0        Sqrt of Final Impact Energy in ev^{1/2}       ' date ])
elseif tval == 3  
  if dval == 2
    semilogy(U/1e3,Efq+1),grid  
    xlabel('Peak voltage   [kV]')
    ax = axis;
    axis([min(U)/1e3,max(U)/1e3,0,ax(4)]);    
  elseif dval == 3
    semilogy(Efl/1e3,Efq+1),grid  
    xlabel('Peak Electric Field   [kV/m]')
    ax = axis;
    axis([min(Efl)/1e3,max(Efl)/1e3,0,ax(4)]);    
  elseif dval == 1
    semilogy(Pow/1e3,Efq+1),grid  
    xlabel('Forward Power   [kW]')
    ax = axis;
    axis([min(Pow)/1e3,max(Pow)/1e3,0,ax(4)]);    
  end
  ylabel(['Ef_{' num2str(N) '} ' ])
  title(['MultiPac 2.0            Final Impact Energy [ev] in log_{10}-scale         ' date ])  
end

if tval == 1 | tval == 3
%  sec1 = secy1(e1,1);
  sec1 = interp1(secy1(1:e1,2),secy1(1:e1,1),1);
  sec2 = secy1(e2,1);
  sec3 = secy1(e3,1);
elseif tval == 2  
%  sec1 = sqrt(secy1(e1,1));
  sec1 = sqrt(interp1(secy1(1:e1,2),secy1(1:e1,1),1));
  sec2 = sqrt(secy1(e2,1));
  sec3 = sqrt(secy1(e3,1));  
end

if dval == 2
  fl1 = min(U)/1e3;
  fl2 = max(U)/1e3;
elseif dval == 3
  fl1 = min(Efl)/1e3;
  fl2 = max(Efl)/1e3;
elseif dval == 1
  fl1 = min(Pow)/1e3;
  fl2 = max(Pow)/1e3;
end

% plot the secondaty yield tresholds and the maxium yield
hold on
plot([fl1,fl2],[sec1,sec1],'-r')
plot([fl1,fl2],[sec2,sec2],'-r')
plot([fl1,fl2],[sec3,sec3],'--r')
hold off

% plot the other secondary yield curves
if exist('secy2')
  load secy2;
  e1 = min(find(secy2(:,2)>=1));      % lower thershold
  e2 = max(find(secy2(:,2)>=1));      % upper thershold
  [val,e3] = max(secy2(:,2));         % maximum secondary yield
  if tval == 1 | tval == 3
    sec1 = interp1(secy2(1:e1,2),secy2(1:e1,1),1);
    sec2 = secy2(e2,1);
    sec3 = secy2(e3,1);
  elseif tval == 2  
    sec1 = sqrt(interp1(secy2(1:e1,2),secy2(1:e1,1),1));
    sec2 = sqrt(secy2(e2,1));
    sec3 = sqrt(secy2(e3,1));  
  end

  % plot the secondaty yield tresholds and the maxium yield
  hold on
  plot([fl1,fl2],[sec1,sec1],'-k')
  plot([fl1,fl2],[sec2,sec2],'-k')
  plot([fl1,fl2],[sec3,sec3],'--k')
  hold off
end  

if exist('secy3')
  load secy3;
  e1 = min(find(secy3(:,2)>=1));      % lower thershold
  e2 = max(find(secy3(:,2)>=1));      % upper thershold
  [val,e3] = max(secy3(:,2));         % maximum secondary yield
  if tval == 1 | tval == 3
    sec1 = interp1(secy3(1:e1,2),secy2(1:e1,1),1);
    sec2 = secy3(e2,1);
    sec3 = secy3(e3,1);
  elseif tval == 2  
    sec1 = sqrt(interp1(secy3(1:e1,2),secy2(1:e1,1),1));
    sec2 = sqrt(secy3(e2,1));
    sec3 = sqrt(secy3(e3,1));  
  end

  % plot the secondaty yield tresholds and the maxium yield
  hold on
  plot([fl1,fl2],[sec1,sec1],'-g')
  plot([fl1,fl2],[sec2,sec2],'-g')
  plot([fl1,fl2],[sec3,sec3],'--g')
  hold off
end  

if exist('secy4')
  load secy4;
  e1 = min(find(secy4(:,2)>=1));      % lower thershold
  e2 = max(find(secy4(:,2)>=1));      % upper thershold
  [val,e3] = max(secy4(:,2));         % maximum secondary yield
  if tval == 1 | tval == 3
    sec1 = interp1(secy4(1:e1,2),secy2(1:e1,1),1);
    sec2 = secy4(e2,1);
    sec3 = secy4(e3,1);
  elseif tval == 2  
    sec1 = sqrt(interp1(secy4(1:e1,2),secy2(1:e1,1),1));
    sec2 = sqrt(secy4(e2,1));
    sec3 = sqrt(secy4(e3,1));  
  end

  % plot the secondaty yield tresholds and the maxium yield
  hold on
  plot([fl1,fl2],[sec1,sec1],'-m')
  plot([fl1,fl2],[sec2,sec2],'-m')
  plot([fl1,fl2],[sec3,sec3],'--m')
  hold off
end  
% --------------------------------------------------------------------
