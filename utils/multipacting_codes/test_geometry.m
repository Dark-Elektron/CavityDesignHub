% Function program ok = test_geometry(s)
% --------------------------------------------------------------------
% To check that the geometry files are defined correctly. So far, only
% a limited number of checks are made.
% --------------------------------------------------------------------
% CALLS TO : error_message.m, clear_window.m
% 17/05/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% 16/06/00 :                 - a bug fixed
% --------------------------------------------------------------------

function ok = test_geometry(s)

if nargin < 1
  s = 1;
end

ok  = 1;
ok1 = exist('geodata.n');
ok2 = exist('geodatal.n');
ok3 = exist('geodataw.n');
ok4 = exist('geodatar.n');

ok5 = check_fieldparam;
if ok5 > 0
  load fieldparam;
  gtype = fieldparam;
else
  gtype = 0;
end

cl = 0;
% ---------------------- cavity and coupler ------------------------
if ok1 > 0 & gtype < 3
  
  if s > 0
    cl = error_message('Testing file geodata.n.');
  end
  
  load geodata.n
  ng = size(geodata,1);
  gr = geodata(4:ng,1);
  gz = geodata(4:ng,2);
  gb = geodata(4:ng,3);
  gs = geodata(4:ng,4);
  nr = length(gr);
  
  [zmin, zind] = min(gz);
  
  tol = min(sqrt((gr(2:nr)-gr(1:nr-1)).^2+(gz(2:nr)-gz(1:nr-1)).^2)) / 1000;
  if tol < eps, tol = 100*eps; end
  
  % check that the first point is at the axis (cavity)
  if gr(1) > tol
    if gtype == 1  
      cl=error_message('Error: For a cavity first point must be at the axis.');
      ok = 0;
    elseif gtype == 0
      cl = error_message(['Warning: For a cavity first point must be at ' ...
			  'the axis.']);
    end    
  end

  % check that the first point is on the left hand side (cavity, coupler)
  if gtype < 3
    if any(zind) ~= 1
      cl = error_message('Error: First point must be on the left hand side.');
      ok = 0;
    end
  end
  
  % check that the first and last point coincide (cavity and coupler)
  dis = sqrt((gr(nr)-gr(1)).^2+(gz(nr)-gz(1)).^2);
  if dis > tol
    cl = error_message('Error: First and last points must coincide.');
    ok = 0;
  end   
  
  % check that the secondary yield corresponds to the boundaries
  SPE = find(gb==0);
  PEC = find(gb==1);
  WIN = find(gb==2);
  PMC = find(gb==3);
  
  sec0 = find(gs==0);
  sec1 = find(gs==1);
  sec2 = find(gs==2);
  sec3 = find(gs==3);
  sec4 = find(gs==4);
  
  if length([SPE(:)',PMC(:)']) ~= length(sec0(:)')
    cl = error_message(['Warning: Non-zero secondary yield for magnetic '...
			'walls or artificial boundaries.']);
%    ok = 0;
  else
    tmp1 = sort([SPE(:)',PMC(:)']);
    tmp2 = sort(sec0(:)');
    if tmp1 ~= tmp2
      cl = error_message(['Warning: Non-zero secondary yield for magnetic ' ...
			  'walls or artificial boundaries.']);
%      ok = 0;
    end
  end
end

% ---------------------- window ------------------------
if ok2*ok3*ok4 > 0
  if s > 0
    cl = error_message('Testing files geodatal.n, geodataw.n and geodatar.n.');
  end

  load geodatal.n
  ngl = size(geodatal,1);
  grl = geodatal(4:ngl,1);
  gzl = geodatal(4:ngl,2);
  gbl = geodatal(4:ngl,3);
  gsl = geodatal(4:ngl,4);
  nrl = length(grl);
  
  load geodataw.n
  ngw = size(geodataw,1);
  grw = geodataw(4:ngw,1);
  gzw = geodataw(4:ngw,2);
  gbw = geodataw(4:ngw,3);
  gsw = geodataw(4:ngw,4);
  nrw = length(grw);
  
  load geodatar.n
  ngr = size(geodatar,1);
  grr = geodatar(4:ngr,1);
  gzr = geodatar(4:ngr,2);
  gbr = geodatar(4:ngr,3);
  gsr = geodatar(4:ngr,4);
  nrr = length(grr);
  
  [zmin,minind] = min(gzl);
  [zmax,maxind] = min(gzr);
  
  gr = [grl;grw;grr];
  gz = [gzl;gzw;gzr];
  nr = length(gr);
  
  tol = min(sqrt((gr(2:nr)-gr(1:nr-1)).^2+(gz(2:nr)-gz(1:nr-1)).^2)) / 1000;
  if tol < eps, tol = 100*eps; end

  % check that the first point is on the left hand side (warm side)
  if any(minind) ~= 1
    cl = error_message(['Error: First point on the warm side must be on ' ...
			'the left hand side.']);
    ok = 0;
  end
  if any(maxind) ~= 1
    cl = error_message(['Error: Last point on the cold side must be on ' ...
			'the right hand side.']);
    ok = 0;
  end
  
  % check that the first and last point coincide (warm side)
  dis = sqrt((grl(nrl)-grl(1)).^2+(gzl(nrl)-gzl(1)).^2);
  if dis > tol
    cl = error_message(['Error: First and last points must coincide ' ...
			'(warm side).']);
    ok = 0;
  end   
  dis = sqrt((grr(nrr)-grr(1)).^2+(gzr(nrr)-gzr(1)).^2);
  if dis > tol
    cl = error_message(['Error: First and last points must coincide ' ...
			'(cold side).']);
    ok = 0;
  end   
  
  % check that the secondary yield corresponds to the boundaries
  SPE = find(gbl==0);
  PEC = find(gbl==1);
  WIN = find(gbl==2);
  PMC = find(gbl==3);
  
  sec0 = find(gsl==0);
  sec1 = find(gsl==1);
  sec2 = find(gsl==2);
  sec3 = find(gsl==3);
  sec4 = find(gsl==4);
  
  if length([SPE(:)',PMC(:)']) ~= length(sec0(:)')
    cl = error_message(['Warning: Non-zero secondary yield for magnetic '...
			'walls or artificial boundaries (warm side).']);
  else
    tmp1 = sort([SPE(:)',PMC(:)']);
    tmp2 = sort(sec0(:)');
    if tmp1 ~= tmp2
      cl = error_message(['Warning: Non-zero secondary yield for magnetic ' ...
			  'walls or artificial boundaries (warm side).']);
    end
  end

  SPE = find(gbr==0);
  PEC = find(gbr==1);
  WIN = find(gbr==2);
  PMC = find(gbr==3);
  
  sec0 = find(gsr==0);
  sec1 = find(gsr==1);
  sec2 = find(gsr==2);
  sec3 = find(gsr==3);
  sec4 = find(gsr==4);
  
  if length([SPE(:)',PMC(:)']) ~= length(sec0(:)')
    cl = error_message(['Warning: Non-zero secondary yield for magnetic '...
			'walls or artificial boundaries (cold side).']);
  else
    tmp1 = sort([SPE(:)',PMC(:)']);
    tmp2 = sort(sec0(:)');
    if tmp1 ~= tmp2
      cl = error_message(['Warning: Non-zero secondary yield for magnetic ' ...
			  'walls or artificial boundaries (cold side).']);
    end
  end

%  WIN = find(gbw==2);
%  if length(WIN) ~= length(gbw)
%    error_message(['Error: Window must be defined as a dielectric surface.']);
%    ok = 0;
%  end
  
end

if ok == 1 & s > 0
  cl = error_message('Geometry files are ok.');
end  
if cl == 1
  clear_window;
end  
% --------------------------------------------------------------------
