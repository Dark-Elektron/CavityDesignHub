% Function program iy0 = mapplot1(map,q,par,y00,bo,trajec)
% -------------------------------------------------------------------------
% Pplots D_N, Ea_N or Ef_N mapping.
% INPUT  map : (1,n) vector giving the place-time distance between the 
%              place of initial emission and final place of electron after 
%              N impacts, n is the number of initial electrons
%        q   : type of parameter map, q == 1, map = D_N, q == 2,
%              map = Ea_N and q == 3, map = Ef_N
%        par : (7,1) vector listing the MP parameters
%        y00 : initial sites
%        bo  : the boundary curve
% OUTPUT iy0 : indexes of the chosen initial electrons referring to matrix 
%              y00
% -------------------------------------------------------------------------
% CALLS TO : none
% 01/02/99 : Marko Ukkola - Rolf Nevanlinna Institute (mapplot.m)
% 07/03/00 : Pasi Ylä-Oijala - modifications (mapplot1.m)
% -------------------------------------------------------------------------

function iy0 = mapplot1(D_N,q,par,y00,bo,trajec)

N = par(5);
if q ~= 1
  if q == 2
    ti = ['Ea_'];
  else
    ti = ['Ef_'];
  end
  q = 0.16e-18;
else
  ti = [' Distance map   d_{' num2str(N) '}'];
end

k = 0;
Prz = zeros(size(y00));
Prz(1,:) = y00(1,:);

ky = 1;
yt = zeros(1,size(y00,1));
yt(1) = y00(1,4);

for i = 1:size(y00,1)
  if y00(i,4) == Prz(1,4)
    k = k + 1;
    Prz(k,:) = y00(i,:);
  end
  
  if y00(i,4) < yt(1)
    yt(2:(ky+1)) = yt(1:ky);
    yt(1) = y00(i,4);
    ky = ky + 1;
  else
    if y00(i,4) > yt(ky)
      yt(ky + 1) = y00(i,4);
      ky = ky + 1;      
    else
      if ((y00(i,4) > yt(1) & y00(i,4) < yt(ky)) ...
	& all(y00(i,4) ~= yt))
	pla = sum(yt(1:ky) < y00(i,4));
	yt((pla + 2):(ky+1)) = yt((pla + 1):(ky));
	yt(pla + 1) = y00(i,4);
	ky = ky + 1;
      end        
    end    
  end  
end
yt = yt(1:ky);
Prz = Prz(1:k,:);

map = zeros(length(yt),size(Prz,1));

for i = 1:size(y00,1)
  iy = find(Prz(:,1) == y00(i,1) & Prz(:,2) == y00(i,2));
  ix = find(yt == y00(i,4));
  map(ix,iy) = D_N(i);
end

if all(all(map == -2))
  map = 100*ones(length(yt),size(Prz,1));
else
  map(find(map == -2)) = ones(size(map(find(map == -2))))*max(max(map))*1.1;
end

subplot(2,1,2)
plot(bo(2,:),bo(1,:),'r',Prz(:,2),Prz(:,1),'or')
vec = 1:3:size(Prz,1);
if vec(length(vec)) ~= size(Prz,1)
  vec = [vec size(Prz,1)];
end
dip = 4e-4;
for i = vec
  text(Prz(i,2)+abs(dip*0.5),Prz(i,1)+dip,num2str(i));
  dip = -dip; 
end
xlabel('z axis [m]')
ylabel('r axis [m]')
title('Initial points')

map(1,1) = map(1,2) + 2*eps;
subplot(2,1,1)
pcolor(1:size(Prz,1),yt*par(1)*360,map)
shading flat
colormap(hot)
matr = colormap;
iso = 1.1*max(max(D_N));
caxis([0 iso]);
colorbar; 
colormap(matr.^(1/3));

ylabel('Initial phase [deg]')
xlabel('Place referring to picture below')
title(['MultiPac 2.0          ' ti '              ' date])

iy0 = [];
if trajec == 1
  [x,y] = ginput; x = x(length(x)); 
  iy0 = zeros(1,length(x));
  for i = 1:length(x)
    tmp = abs(yt*par(1)*360 - y(i));
    y(i) = min(find(tmp == min(tmp)));
    xd = [max([round(x(i)-2) 1]) min([round(x(i)+2) size(map,2)])];
    yd = [max([round(y(i)-2) 1]) min([round(y(i)+2) size(map,1)])];
    smap = map(yd(1):yd(2),xd(1):xd(2));
    mi = min(min(smap));
    [J,I] = find(smap == mi);
    indx = xd(1) + I(1) - 1;
    indy = yd(1) + J(1) - 1;
    tmp = find(D_N == mi);
    if map(indy,indx) == mi
      iy0(i) = find(y00(:,1) == Prz(indx,1) & ...
               y00(:,2) == Prz(indx,2) & ...
 	       y00(:,4) == yt(indy));
    end
  end  
end  
% -------------------------------------------------------------------------
