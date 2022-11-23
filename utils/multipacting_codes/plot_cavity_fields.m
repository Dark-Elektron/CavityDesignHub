% ----------------------------------------------------------------------
% Plots the fields computed by the field solver. For a cavity.
%
% INPUT  s : 0 - pcolor plots
%            1 - arrow plots
%  gridcons : grid constant, [mm]
% ----------------------------------------------------------------------
% CALLS TO : arrow.m, error_message.m, clear_window.m
% 17/04/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
% ----------------------------------------------------------------------

function plot_cavity_fields(s,gridcons,ss)

if nargin < 1
  s = 0;
end
if nargin < 3
  ss = 1;
end

cl = 0;
load geodata.n
if nargin < 2
  gridcons = geodata(1,1);
end  
n  = length(geodata(:,1));
gr = geodata(4:n,1);
gz = geodata(4:n,2);

% load the field values
load fields Er Ez H I1 I2 zz rr z r
  
% plot the fields
if s == 0                                   % a pcolor plot
  Er = reshape(Er,length(I2),length(I1));
  Ez = reshape(Ez,length(I2),length(I1));
  EE = sqrt(abs(Er).^2+abs(Ez).^2);
  
  % Calculate accerating field
  L_half_cell = 0.187; % length of cavity half cell
  E_axis = EE(rr==0); %& zz>=-L_half_cell & zz <= L_half_cell
  z_slice = zz(1,:);
  z_active = z_slice(z_slice>=-L_half_cell & z_slice <= L_half_cell);
  Vacc = trapz(z_slice, E_axis);
  Eacc = Vacc/(2*L_half_cell);
  
  Epk_Eacc = max(max(EE))/Eacc
  
  subplot(2,1,1)
  pcolor(zz,rr,EE),colorbar
  shading interp, colormap jet
  hold on
  plot(gz,gr,'-b')
  hold off
  xlabel('z axis [m]')
  ylabel('r axis [m]')
  title(['MultiPac 2.0          Electric field   abs (E)   [V/m]      ' date ])
  axis equal  

  HH = reshape(H,length(I2),length(I1));
  subplot(2,1,2)
  pcolor(zz,rr,4e-7*pi*HH); colorbar
  shading interp, colormap jet
  hold on
  plot(gz,gr,'-b')
  hold off
  title('Magnetic field     B_\phi  [TESLA]')
  xlabel('z axis [m]')
  ylabel('r axis [m]')
  axis equal

else                                   % an arrow plot
  F = [Ez'; Er'];
  A = [z';r'];
  
  subplot(1,1,1)
%   arrow(F,A,gridcons,'f','s','k')
  quiver(z, r, Ez, Er);
  hold on
  plot(gz,gr,'-b')
%  xlabel('z axis [m]')
  ylabel('r axis [m]')
  hold off
  title(['MultiPac 2.0            Electric field        [V/m]         ' date ])
  axis equal
  
  if ss > 0
    cl = error_message(['Electric field plotted. To see the magnetic ' ...
			'field choose Pcolor.']);
  end
end

if cl == 1
  clear_window;
end  
% ----------------------------------------------------------------------
