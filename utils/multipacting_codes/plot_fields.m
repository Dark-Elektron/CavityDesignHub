% Function program plot_fields(s,gridcons,ERr,EZr,BPr,ERi,EZi,BPi,r,z)
% ---------------------------------------------------------------------
% Plot the mixed electromagnetic fields.
% INPUT  s : 0 - pcolor
%            1 - arrow
%        gridcons : grid constant for plotting arrow plots
% ---------------------------------------------------------------------
% CALLS TO : arrow.m, error_message.m
% 18/05/00 : Pasi Yla-Oijala - modifications
% 17/10/00 :                 - coefficient files removed
% ---------------------------------------------------------------------

function plot_fields(s,gridcons,ERr,EZr,BPr,ERi,EZi,BPi,r,z)

load fieldparam;
gtype = fieldparam(1);

if gtype <= 2
  load geodata.n
  n  = length(geodata(:,1));
  gr = geodata(4:n,1);
  gz = geodata(4:n,2);
elseif gtype == 3
  load geodatal.n
  load geodataw.n
  load geodatar.n

  nl  = length(geodatal(:,1));
  grl = geodatal(4:nl,1);
  gzl = geodatal(4:nl,2);
  nw  = length(geodataw(:,1));
  grw = geodatar(4:nw,1);
  gzw = geodatar(4:nw,2);
  nr  = length(geodatar(:,1));
  grr = geodatar(4:nr,1);
  gzr = geodatar(4:nr,2);
end    

if s == 0
  EEr = sqrt(abs(ERr).^2+abs(EZr).^2);
  subplot(2,2,1)
  pcolor(z,r,EEr)
  colormap hot
  hold on
  if gtype <= 2
    plot(gz,gr,'-b')
  else
    plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  end
  hold off
  axis equal
  colorbar,shading interp  
  title('MultiPac 2.0        |Real(E)|   [V/M]')
  ylabel('r axis [m]')
  xlabel('z axis [m]')
      
  EEi = sqrt(abs(ERi).^2+abs(EZi).^2);
  subplot(2,2,2)
  pcolor(z,r,EEi)
  colormap hot
  hold on
  if gtype <= 2
    plot(gz,gr,'-b')
  else
    plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  end
  hold off
  axis equal
  colorbar,shading interp  
  title(['|Imag(E)|   [V/M]         ' date ])
  ylabel('r axis [m]')
  xlabel('z axis [m]')
      
  subplot(2,2,3)
  pcolor(z,r,BPr)
  colormap hot
  colorbar,shading interp  
  hold on
  if gtype <= 2
    plot(gz,gr,'-b')
  else
    plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  end
  hold off
  axis equal
  title('Real(B_\phi) [TESLA]')
  ylabel('r axis [m]')
  xlabel('z axis [m]')

  subplot(2,2,4)
  pcolor(z,r,BPi)
  colormap hot
  colorbar,shading interp  
  hold on
  if gtype <= 2
    plot(gz,gr,'-b')
  else
    plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  end
  hold off
  axis equal
  title('Imag(B_\phi) [TESLA]')
  ylabel('r axis [m]')
  xlabel('z axis [m]')

else
  Fr = [EZr(:),ERr(:)].';
  Fi = [EZi(:),ERi(:)].';
  tmp = points3(r(:)',0,z(:)');
  A = tmp([3,1],:);

  subplot(2,1,1)
  if max(max(abs(Fr))) > 0
    arrow(Fr,A,gridcons,'f','s','k')
    hold on
  end
  if gtype <= 2
    plot(gz,gr,'-b')
  else
    plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  end
  xlabel('z axis [m]')
  ylabel('r axis [m]')
  hold off
  axis equal
  title(['MultiPac 2.0            Real(E)    [V/m]            ' date ])
  axis equal
      
  subplot(2,1,2)
  if max(max(abs(Fi))) > 0
    arrow(Fi,A,gridcons,'f','s','k')    
    hold on
  end
  if gtype <= 2
    plot(gz,gr,'-b')
  else
    plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  end
  xlabel('z axis [m]')
  ylabel('r axis [m]')
  hold off
  axis equal
  title('Imag(E)    [V/m]')
  axis equal

  if s > 0
    error_message(['Electric field plotted. To see the magnetic field ' ...
		   'choose Pcolor.']);
  end
end
% ---------------------------------------------------------------------
