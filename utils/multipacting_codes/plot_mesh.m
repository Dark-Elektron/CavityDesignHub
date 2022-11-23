% Function program plot_mesh(s)
% --------------------------------------------------------------------
% Plots the mesh.
%
% --------------------------------------------------------------------
% CALLS TO : error_message.m, clear_window.m
% xx/yy/00 : Seppo J‰rvemp‰‰
% 04/04/00 : Pasi Yl‰-Oijala - m-file
% --------------------------------------------------------------------

function plot_mesh(s)


if nargin < 1
  s = 1;
end

cl = 0;
ok1 = exist('mesh.mat');
ok2 = exist('fieldparam');
if ok1 == 0
  cl = error_message(['The mesh does not exists. Choose Mesh Generator ' ...
		      'in menu Run.']);
elseif ok2 == 0
  cl = error_message('Parameter file fieldparam does not exist.');
else  

  load fieldparam;
  gtype = fieldparam(1);
  if gtype == 1        
    if s > 0
      cl = error_message('Plotting the mesh.');
%      error_message('                  ');
    end
  else
    if s > 0
      cl = error_message('Plotting the mesh. Blue area for streching.');
%      error_message('                  ');
    end
  end
  
  Window = figure(3);
  clf;
  set(Window,'name','MULTIPAC - Input Window II');

  % plots 2-d mesh in mesh.mat, which includes coord and etopol -arrays
  load mesh

  xmin = min(coord(:,1));
  xmax = max(coord(:,1));
  ymin = min(coord(:,2));
  ymax = max(coord(:,2));
  plot([xmin xmax xmax xmin], [ymin ymin ymax ymax], 'k.')
  hold on

  [koko, pois] = size(etopol);
  ala = 0.0;
  X = zeros(4, koko); Y = zeros(4, koko);
  C = zeros(1, koko);
  for i = 1:koko
    n1 = etopol(i,1); x1 = coord(n1, 1); y1 = coord(n1, 2);
    n2 = etopol(i,2); x2 = coord(n2, 1); y2 = coord(n2, 2);
    n3 = etopol(i,3); x3 = coord(n3, 1); y3 = coord(n3, 2);
    X(:,i) = [x1 x2 x3 x1]';
    Y(:,i) = [y1 y2 y3 y1]';
  
    osa = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
    ala = ala + osa;
    plot([x1 x2 x3 x1], [y1 y2 y3 y1], 'k');
  end

  I = find(alue == 0);
  fill(X(:,I), Y(:,I),'b');
  I = find(tyyppi == 1);
  fill(X(:,I), Y(:,I),'r');
  I = find(tyyppi == 2);
  fill(X(:,I), Y(:,I),'g');

  I = find(edges(:,3) > 0); % no bouncing
  for i = 1:length(I)
    i1 = edges(I(i),1); i2 = edges(I(i),2);
    plot([coord(i1,1), coord(i2,1)], [coord(i1,2), coord(i2,2)], 'r')
  end
  I = find(boundary(1,4) == 3);
  for i = 1:length(I)
    plot(coord(I(i),1), coord(I(i),2), 'b*')
  end
  I = find(boundary(1,4) == 0);
  plot(coord(I,1), coord(I,2), 'w*')

  ala = ala/2;
  title(['MultiPac 2.0                    Mesh                 ' date ])
  hold off
  xlabel('z axis [m]')
  ylabel('r axis [m]')
  axis('image')
  colormap(jet)
end  

if cl == 1
  clear_window;
end  
% --------------------------------------------------------------------

