% Function program plot_trajectory(side,s)
% ---------------------------------------------------------------------------
% Plot an electron trajectory.
% INPUT  side : 'l'eft or 'r'ight for a window, 'c' otherwise
%
% ---------------------------------------------------------------------------
% CALLS TO : check_inputs_win, check_inputs, error_message.m, clear_window.m
%            check_geotype.m
% 09/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% 13/12/00 :                 - check_geotype.m added
% ---------------------------------------------------------------------------

function plot_trajectory(side,s)

if nargin < 1
  side = 'c';
end
if nargin < 2
  s = 1;
end

ok3 = 0;
ok3 = check_fieldparam;
if ok3 > 0
  load fieldparam
  gtype = fieldparam(1);
  if side == 'l' | side == 'r'
    if gtype < 3
      if s > 0
	cl = error_message(['The geometry is not a window. Choose Cavity '...
			    'or Coupler.']);
      end
      ok3 = 0;
    end
  end
  if side == 'c'
    if gtype == 3
      if s > 0
	cl = error_message('The geometry is a window. Choose Window.');
      end
      ok3 = 0;
    end
  end  
end
%ok3 = check_geotype(side,0);
if ok3 == 0
  return;
end

cl = 0;
if side == 'l' | side == 'r'
  ok1 = check_inputs_win;
else  
  ok1 = check_inputs;
end  
ok2 = exist('elecpath');

if (ok1*ok2 == 0)
  if ok2 == 0
    cl=error_message(['Electron trajectory is missing. Choose Trajectory ' ... 
		      'in menu Run.']);
  end
else  

  load fieldparam
  gtype = fieldparam(1);
  if gtype < 3
    if s > 0
      cl = error_message('Plotting an electron trajectory.');
    end
  elseif gtype == 3
    if s > 0 & (side == 'l' | side == 'r')
      cl = error_message('Plotting an electron trajectory in window.');
    end
    if side == 'c'
      cl = error_message(['To plot an electron trajectory in window, ' ...
			  'choose warm or cold side.']);
      return;
    end  
  end    
  
  Window = figure(7);
  clf;
  set(Window,'name','MULTIPAC - Output Window III');

  if side == 'l' | side == 'r'
    load geodatal.n
    ng  = length(geodatal(:,1));
    bol = geodatal(4:ng,1:2)';
    load geodatar.n
    ng  = length(geodatar(:,1));
    bor = geodatar(4:ng,1:2)';
    bo  = [bol,bor];
    load geodataw.n
    ng  = length(geodataw(:,1));
    wr  = geodataw(4:ng,1);
    wz  = geodataw(4:ng,2);
  else
    load geodata.n
    ng  = length(geodata(:,1));
    bo  = geodata(4:ng,1:2)';
    wr  = [];
    wz  = [];
  end
  load param
  load elecpath
  par = param;

  n = size(elecpath,1);
  if (n == 1)
    pat = [];
    cl = error_message(['No electron emission. Please, define a new ' ...
			'initial point.']);
  else
    figure(7)
    pat = elecpath(2:n,[1 3 4 6 7 8]);
    N = par(5);
    hit = find(pat(:,6) ~= 0);
    hit = hit(2:2:length(hit));
    speed = sqrt(pat(hit,3).^2 + pat(hit,4).^2);
    c = 2.9979e8; 
    M = 9.1093879e-31; 
    q = 1.6021773e-19;
    energy = (1./sqrt(1.0-(speed.^2)./(c.^2)) - 1) .* M .* c.^2 ./ q;
    avegene = num2str(mean(energy));
    finaene = num2str(energy(length(energy)));
    maxiene = num2str(max(energy));

    subplot(3,1,1)
    if side == 'l' | side == 'r'
      plot(bol(2,:),bol(1,:),'-b',bor(2,:),bor(1,:),'-b')
      hold on
      fill(wz,wr,'y');
      hold off
    else
      plot(bo(2,:),bo(1,:),'-b')
    end
    hold on
    plot(pat(:,2),pat(:,1),'-r')
    hold off
    title(['MultiPac 2.1       Electron Trajectory,   N = ' num2str(N) ...
	   ',     ' date])
    dt = abs(min(pat(:,5)) - max(pat(:,5))) * par(1);
    xlabel(['z-axis [m],  flight time ' num2str(dt) ' periods'])
    ylabel('r-axis [m]')
    grid
    axis equal

    subplot(3,1,2)
    if side == 'l' | side == 'r'
      plot(bol(2,:),bol(1,:),'-b',bor(2,:),bor(1,:),'-b')      
      hold on
      fill(wz,wr,'y');
      hold off
    else
      plot(bo(2,:),bo(1,:),'-b')
    end
    hold on
    plot(pat(:,2),pat(:,1),'-r',pat(hit,2),pat(hit,1),'ro')
    hold off
    min1 = 0.9*min(pat(:,1))-eps;
    max1 = 1.1*max(pat(:,1))+eps;
    if min(pat(:,2)) < 0
      min2 = 1.1*min(pat(:,2))-eps;
    else
      min2 = 0.9*min(pat(:,2))-eps;
    end      
    if max(pat(:,2)) < 0
      max2 = 0.9*max(pat(:,2))+eps;
    else
      max2 = 1.1*max(pat(:,2))+eps;
    end      
    axis([min2,max2,min1,max1])
    grid
    xlabel('z-axis [m]')
    ylabel('r-axis [m]')

    subplot(3,1,3)
    plot(pat(:,5)*par(1),pat(:,1),'r',pat(hit,5)*par(1),pat(hit,1),'ro'),
    midax = axis;
    midax(3) = max(min(bo(1,:)),min(pat(:,1) - 1e-3));
    midax(4) = min(max(bo(1,:)),max(pat(:,1) + 1e-3));
    axis(midax)
   xlabel(['time in [1/f], average energy ' avegene ', final energy ' finaene])
    ylabel('r-axis [m]')

%    subplot(3,1,3)
%    plot(pat(:,5)*par(1),pat(:,2),'r',pat(hit,5)*par(1),pat(hit,2),'ro')
%    botax = axis;
%    botax(3) = max(min(bo(2,:)),min(pat(:,2) - 1e-3));
%    botax(4) = min(max(bo(2,:)),max(pat(:,2) + 1e-3));
%    axis(botax)
%    ylabel('z-axis')
%    xlabel('time in [1/f]')
  end    
end  

if cl == 1
  clear_window;
end
% --------------------------------------------------------------------


