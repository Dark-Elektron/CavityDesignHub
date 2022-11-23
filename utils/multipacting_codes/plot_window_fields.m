% Function program plot_window_fields(s,gridcons,ss)
% ----------------------------------------------------------------------
% Plots the fields computed by the field solver. For a window.
%
% INPUT  s : 0 - pcolor plots
%            1 - arrow plots
%  gridcons : grid constant, [mm]
% ----------------------------------------------------------------------
% CALLS TO : arrow.m, error_message.m, clear_window.m
% 17/04/00 : Pasi Ylä-Oijala - Rolf Nevalinna Institute
% ----------------------------------------------------------------------

function plot_window_fields(s,gridcons,ss)

if nargin < 1
  s = 0;
end
if nargin < 3
  ss = 1;
end

cl = 0;

% load the geometry
load geodatal.n
nl  = length(geodatal(:,1));
grl = geodatal(4:nl,1);
gzl = geodatal(4:nl,2);
load geodataw.n
nw  = length(geodataw(:,1));
grw = geodataw(4:nw,1);
gzw = geodataw(4:nw,2);
load geodatar.n
nr  = length(geodatar(:,1));
grr = geodatar(4:nr,1);
gzr = geodatar(4:nr,2);

if nargin < 2
  gridcons = geodatal(1,1);
end  

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
  
% plot the fields
if s == 0                                   % a pcolor plot
  Er = reshape(Er1,length(I21),length(I11));
  Ez = reshape(Ez1,length(I21),length(I11));
  EE = sqrt(abs(Er).^2+abs(Ez).^2);
  
  subplot(2,2,1)
  pcolor(zz1,rr1,EE),colorbar
  shading interp, colormap hot
  hold on
  plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  hold off
  xlabel('z axis [m]             E-walls')
  ylabel('r axis [m]')
  title('MultiPac 2.0      |E|   [V/m]')
  axis equal
  
  Er = reshape(Er2,length(I22),length(I12));
  Ez = reshape(Ez2,length(I22),length(I12));
  EE = sqrt(abs(Er).^2+abs(Ez).^2);
  
  subplot(2,2,3)
  pcolor(zz2,rr2,EE),colorbar
  shading interp, colormap hot
  hold on
  plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  hold off
  xlabel('z axis [m]            H-walls')
  ylabel('r axis [m]')
  title('Electric field    |E|   [V/m]')
  axis equal
  
  HH = reshape(H1,length(I21),length(I11));
  subplot(2,2,2)
  pcolor(zz1,rr1,4e-7*pi*HH); colorbar
  shading interp, colormap hot
  hold on
  plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  hold off
  title(['B_\phi   [TESLA]       ' date ])
  xlabel('z axis [m]             E-walls')
  ylabel('r axis [m]')
  axis equal

  HH = reshape(H2,length(I22),length(I12));
  subplot(2,2,4)
  pcolor(zz2,rr2,4e-7*pi*HH); colorbar
  shading interp, colormap hot
  hold on
  plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  hold off
  title('Magnetic field     B_\phi   [TESLA]')
  xlabel('z axis [m]            H-walls')
  ylabel('r axis [m]')
  axis equal

else                                   % an arrow plot
  F = [Ez1'; Er1'];
  A = [z1';r1'];

  subplot(2,1,1)
  fill(gzw,grw,'y')
  hold on
  arrow(F,A,gridcons,'f','s','k')
  plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
  hold off
  xlabel('z axis [m]            E-walls')
  ylabel('r axis [m]')
  title(['MultiPac 2.0             Electric field        [V/m]     ' date ])
  axis equal

  F = [Ez2'; Er2'];
  A = [z2';r2'];

  subplot(2,1,2)
  fill(gzw,grw,'y')
  hold on
  arrow(F,A,gridcons,'f','s','k')
  plot(gzl,grl,'-b',gzw,grw,'-b',gzr,grr,'-b')
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

