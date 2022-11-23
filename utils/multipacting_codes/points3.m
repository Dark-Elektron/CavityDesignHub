% FUNCTION program X = POINTS3(x,y,z)
% -----------------------------------------------------------------------------
% Program to generate a field-point-matrix X in R^3 for given coordinate- 
% vectors x,y and z. Each column of X is a point in R^3 as follows:
%              X(:,h) = [x(i),y(j),z(k)]', h=1,...,nx*ny*nz.
% The indexing system is as follows:
%              leading index k runs from 1 to nz,
%              for every k, the secondary index j runs from 1 to ny and
%              for every j, the index i runs from 1 to nx.
% The vectors x, y and z can be given in rectangular, spherical (x=r,y=theta,
% z=phi) or cylindrical (x=r,y=phi,z=z) coordinates.
% INPUT   x,y,z : (1,nx), (1,ny) and (1,nz)-vectors specifying the x-,y- and
%                  z-coordinates, [m] or [rad]
% OUTPUT    X   : (3,nx*ny*nz)-matrix giving the field points, [m] or [rad]
% ---------------------------------------------------------------------------
% CALLS TO : None
% 08/10/91 : Pasi Yla-Oijala: Rolf Nevanlinna Institute
% ---------------------------------------------------------------------------

function X = points3(x,y,z)

nx = length(x) ;
ny = length(y) ;
nz = length(z) ;
X = zeros(3,nx*ny*nz) ;

if nz == 1
  zx = ones(1,nx)*z ;
  for j=1:ny
    yx = ones(1,nx)*y(j) ;
    X(:,(j-1)*nx+1:j*nx) = [x;yx;zx] ;
  end
elseif ny == 1
  yx = ones(1,nx)*y ;
  for k = 1:nz
    zx = ones(1,nx)*z(k) ;
    X(:,(k-1)*nx+1:k*nx) = [x;yx;zx] ;
  end
elseif nx == 1
  xx = ones(1,ny)*x ;
  for k=1:nz
    zy = ones(1,ny)*z(k) ;
    X(:,(k-1)*ny+1:k*ny) = [xx;y;zy] ;
  end
else
  ind = 1 ;
  for k = 1:nz
    zx = ones(1,nx)*z(k) ;
    for j = 1:ny
      yx = ones(1,nx)*y(j) ;
      X(:,nx*(ind-1)+1:nx*ind) = [x;yx;zx] ;
      ind = ind + 1 ;
    end
  end
end
% -----------------------------------------------------------------------
