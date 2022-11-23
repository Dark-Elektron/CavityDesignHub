% Function program plot_geometry.m
% -------------------------------------------------------------------------
% Plots the geometry.
%   - Artificial boundaries are black.
%   - Metallic boundaries are red.
%   - Dielectric surfaces are green.
%   - Magnetic walls are blue.
%   - Window is denoted by yellow.
% -------------------------------------------------------------------------
% CALLS TO : check_geodata.m, error_message.m
% 10/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% -------------------------------------------------------------------------

function plot_geometry

ok = check_geodata;

if ok > 0
  load fieldparam
  gtype = fieldparam(1);
  
  cl = error_message('Plotting the geometry.');

  Window = figure(3);
  clf;
  set(Window,'name','MULTIPAC - Input Window II');

  subplot(1,1,1)

  if gtype <= 2
    load geodata.n
    n  = length(geodata(:,1));
    gr = geodata(4:n,1);
    gz = geodata(4:n,2);

    SPE = find(geodata(4:n,3)==0);
    PEC = find(geodata(4:n,3)==1);
    WIN = find(geodata(4:n,3)==2);
    PMC = find(geodata(4:n,3)==3);

    ns  = length(SPE); 
    if ns > 0
      if SPE(ns) == n-3, ns = ns-1; end
    end
    np  = length(PEC); 
    if np > 0
      if PEC(np) == n-3, np = np-1; end
    end
    nm  = length(PMC); 
    if nm > 0
      if PMC(nm) == n-3, nm = nm-1; end
    end
    
    for j=1:ns
      plot([gz(SPE(j)),gz(SPE(j)+1)],[gr(SPE(j)),gr(SPE(j)+1)],'-k'),hold on
    end
    for j=1:np
      plot([gz(PEC(j)),gz(PEC(j)+1)],[gr(PEC(j)),gr(PEC(j)+1)],'-r'),hold on
    end
    for j=1:nm
      plot([gz(PMC(j)),gz(PMC(j)+1)],[gr(PMC(j)),gr(PMC(j)+1)],'-b'),hold on
    end
    hold off
    
  else
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

    SPE = find(geodatal(4:nl,3)==0);
    PEC = find(geodatal(4:nl,3)==1);
    WIN = find(geodatal(4:nl,3)==2);
    PMC = find(geodatal(4:nl,3)==3);

    ns  = length(SPE); 
    if ns > 0
      if SPE(ns) == nl-3, ns = ns-1; end
    end
    np  = length(PEC); 
    if np > 0
      if PEC(np) == nl-3, np = np-1; end
    end
    nw  = length(WIN); 
    if nw > 0
      if WIN(nw) == nl-3, nw = nw-1; end
    end
    nm  = length(PMC); 
    if nm > 0
      if PMC(nm) == nl-3, nm = nm-1; end
    end

    fill(gzw,grw,'y')
    hold on
    plot(gzw,grw,'-g')

    for j=1:ns
      plot([gzl(SPE(j)),gzl(SPE(j)+1)],[grl(SPE(j)),grl(SPE(j)+1)],'-k')
    end
    for j=1:np
      plot([gzl(PEC(j)),gzl(PEC(j)+1)],[grl(PEC(j)),grl(PEC(j)+1)],'-r')
    end
    for j=1:nw
      plot([gzl(WIN(j)),gzl(WIN(j)+1)],[grl(WIN(j)),grl(WIN(j)+1)],'-g')
    end
    for j=1:nm
      plot([gzl(PMC(j)),gzl(PMC(j)+1)],[grl(PMC(j)),grl(PMC(j)+1)],'-b')
    end

    SPE = find(geodatar(4:nr,3)==0);
    PEC = find(geodatar(4:nr,3)==1);
    WIN = find(geodatar(4:nr,3)==2);
    PMC = find(geodatar(4:nr,3)==3);

    ns  = length(SPE); 
    if ns > 0
      if SPE(ns) == nr-3, ns = ns-1; end
    end
    np  = length(PEC); 
    if np > 0
      if PEC(np) == nr-3, np = np-1; end
    end
    nw  = length(WIN); 
    if nw > 0
      if WIN(nw) == nr-3, nw = nw-1; end
    end
    nm  = length(PMC); 
    if nm > 0
      if PMC(nm) == nr-3, nm = nm-1; end
    end

    for j=1:ns
      plot([gzr(SPE(j)),gzr(SPE(j)+1)],[grr(SPE(j)),grr(SPE(j)+1)],'-k')
    end
    for j=1:np
      plot([gzr(PEC(j)),gzr(PEC(j)+1)],[grr(PEC(j)),grr(PEC(j)+1)],'-r')
    end
    for j=1:nw
      plot([gzr(WIN(j)),gzr(WIN(j)+1)],[grr(WIN(j)),grr(WIN(j)+1)],'-g')
    end
    for j=1:nm
      plot([gzr(PMC(j)),gzr(PMC(j)+1)],[grr(PMC(j)),grr(PMC(j)+1)],'-b')
    end
    
    hold off
  end 
  
  axis equal
  grid
  xlabel('z axis [m]')
  ylabel('r axis [m]')
  title(['MultiPac 2.0                    Geometry                ' date ])
end  
% ---------------------------------------------------------------------
