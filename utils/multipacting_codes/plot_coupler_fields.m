% Function program plot_coupler_fields(s,gridcons,ss)
% ----------------------------------------------------------------------
% Plots the fields computed by the field solver. For a coupler.
%
% INPUT  s : 0 - pcolor plots
%            1 - arrow plots
%  gridcons : grid constant, [mm]
% ----------------------------------------------------------------------
% CALLS TO : arrow.m, error_message.m, clear_window.m
% 17/04/00 : Pasi Ylä-Oijala - Role Nevanlinna Institute
% ----------------------------------------------------------------------

function plot_coupler_fields(s,gridcons,ss)

if nargin < 3
  ss = 1;
end

load geodata.n
gridcons = geodata(1,1);
  
% load the geometry
load geodata.n
n  = length(geodata(:,1));
gr = geodata(4:n,1);
gz = geodata(4:n,2);
  
% load the field values
load fields1
Er1 = Er; Ez1 = Ez; H1 = H; 
I11 = I1; I21 = I2; 
zz1 = zz; rr1 = rr;
z1  = z;  r1  = r;

load fields2
Er2 = Er; Ez2 = Ez; H2 = H;
I12 = I1; I22 = I2; 
zz2 = zz; rr2 = rr;
z2  = z;  r2  = r;

cl = 0;
% plot the fields
if s == 0                                   % a pcolor plot
  Er = reshape(Er1,length(I21),length(I11));
  Ez = reshape(Ez1,length(I21),length(I11));
  EE = sqrt(abs(Er).^2+abs(Ez).^2);
  
  subplot(2,2,1)
  pcolor(zz1,rr1,EE),colorbar
  shading interp, colormap hot
  hold on
  plot(gz,gr,'-b')
  hold off
  axis equal
  xlabel('z axis [m]             E-walls')
  ylabel('r axis [m]')
  title('MultiPac 2.0      |E|   [V/m]')
  
  Er = reshape(Er2,length(I22),length(I12));
  Ez = reshape(Ez2,length(I22),length(I12));
  EE = sqrt(abs(Er).^2+abs(Ez).^2);
  
  subplot(2,2,3)
  pcolor(zz2,rr2,EE),colorbar
  shading interp, colormap hot
  hold on
  plot(gz,gr,'-b')
  hold off
  axis equal
  xlabel('z axis [m]            H-walls')
  ylabel('r axis [m]')
  title('Electric field    |E|   [V/m]')
  
  HH = reshape(H1,length(I21),length(I11));
  subplot(2,2,2)
  pcolor(zz1,rr1,4e-7*pi*HH); colorbar
  shading interp, colormap hot
  hold on
  plot(gz,gr,'-b')
  hold off
  axis equal
  title(['B_\phi   [TESLA]       ' date ])
  xlabel('z axis [m]             E-walls')
  ylabel('r axis [m]')

  HH = reshape(H2,length(I22),length(I12));
  subplot(2,2,4)
  pcolor(zz2,rr2,4e-7*pi*HH); colorbar
  shading interp, colormap hot
  hold on
  plot(gz,gr,'-b')
  hold off
  axis equal  
  title('Magnetic field     B_\phi   [TESLA]')
  xlabel('z axis [m]            H-walls')
  ylabel('r axis [m]')

else                                   % an arrow plot

  F = [Ez1'; Er1'];
  A = [z1';r1'];

  subplot(2,1,1)
  arrow(F,A,gridcons,'f','s','k')
  hold on
  plot(gz,gr,'-b')
  hold off
  xlabel('z axis [m]            E-walls')
  ylabel('r axis [m]')
  title(['MultiPac 2.0             Electric field        [V/m]     ' date ])
  axis equal

  F = [Ez2'; Er2'];
  A = [z2';r2'];

  subplot(2,1,2)
  arrow(F,A,gridcons,'f','s','k')
  hold on
  plot(gz,gr,'-b')
  hold off
  xlabel('z axis [m]             H-walls')
  ylabel('r axis [m]')
  title('Electric field        [V/m]')
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
