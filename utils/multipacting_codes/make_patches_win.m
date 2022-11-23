% Function program [nodes,edges,patches] = make_patches_win(geodata,nodes1,
% edges1,releps)
% -------------------------------------------------------------------------
% Converts geodata-formt into (nodes,edges,pacthes)-format. For a window.
% 
% -------------------------------------------------------------------------
% CALLS TO : error_message.m, clear_window.m
% 10/05/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
% 15/06/00 :                 - bugs fixed
% -------------------------------------------------------------------------

function [nodes,edges,patches] = make_patches_win(geodata,nodes1,edges1,releps)

cl = 0;
if releps == 1
  cl = error_message('Relative pervittyvity of the window is 1?');
end

n  = length(geodata(:,1));
gr = geodata(4:n,1);
gz = geodata(4:n,2);
gn = geodata(4:n,3);
n  = length(gr);

nodes   = [];
edges   = [];
patches = [];

nd  = length(nodes1);
tol = min(sqrt((nodes1(2:nd,1)-nodes1(1:nd-1,1)).^2 + ...
	       (nodes1(2:nd,2)-nodes1(1:nd-1,2)).^2)) / 1e4;

% add new nodes
for j=1:n
  dis = min(sqrt((nodes1(:,1)-gz(j)).^2 + (nodes1(:,2)-gr(j)).^2));
  if dis > tol
    nodes = [nodes;[gz,gr]];
    cl = error_message('Warning: Window has new nodes.');
  end
end

% find the edges
n   = length(gr);
nod = [gz,gr];
edg = [(1:n-1)',[2:n-1,1]',gn(1:n-1)];

% nodes due to the edges in the grid
z1  = nodes1(edges1(:,1),1);
r1  = nodes1(edges1(:,1),2);
z2  = nodes1(edges1(:,2),1);
r2  = nodes1(edges1(:,2),2);

zw1 = nod(edg(:,1),1);
rw1 = nod(edg(:,1),2);
zw2 = nod(edg(:,2),1);
rw2 = nod(edg(:,2),2);

% patch for the window
new = max(max(edges1))+1;
for j=1:length(edg)
  [dis1,ind1] = min(sqrt((z1-zw1(j)).^2 + (r1-rw1(j)).^2) + ...
	            sqrt((z2-zw2(j)).^2 + (r2-rw2(j)).^2));  
  [dis2,ind2] = min(sqrt((z1-zw2(j)).^2 + (r1-rw2(j)).^2) + ...
	            sqrt((z2-zw1(j)).^2 + (r2-rw1(j)).^2));  
%  [dis1,ind1,dis2,ind2,j] 
  [dis,jj] = min([dis1,dis2]);
  tmp = [ind1,ind2];
  ind = tmp(jj);

  if dis > tol
    newedge = edg(j,:);
    newnode = nod(newedge(1:2),:);
    [dis1,ind1] = min(sqrt((nodes1(:,1)-newnode(1,1)).^2 + ...
			   (nodes1(:,2)-newnode(1,2)).^2));
    [dis2,ind2] = min(sqrt((nodes1(:,1)-newnode(2,1)).^2 + ...
			   (nodes1(:,2)-newnode(2,2)).^2));
    edges   = [edges;ind1,ind2,2];
    patches = [patches,new];
    new = new+1;
  else
    patches = [patches,ind];
  end
end

patches = [patches,releps,1,1];

if cl == 1
  clear_window;
end  
% ------------------------------------------------------------------------

