% Function program plot_FEM_fields(s,gridcons,ss)
% --------------------------------------------------------------------------
% Plots the fields computed by the field solver.
%
% INPUT  ptype : 0 - pcolor plots, 1 - arrow plots
%  gridcons : grid constant, [mm]
% ---------------------------------------------------------------------------
% CALLS TO : check_geodata.m, check_fields.m, error_message.m, clear_window.m
%            plot_cavity_fields.m, plot_coupler_fields.m, plot_window_fields.m
% 17/04/00 : Pasi Ylä-Oijala - RNI
% ---------------------------------------------------------------------------

function plot_FEM_fields(ptype,gridcons,s)

if nargin < 3
  s = 1;
end

cl  = 0;
ok1 = check_geodata;
ok2 = 0;
if ok1 > 0
  load fieldparam
  gtype = fieldparam(1);
  ok2 = check_fields(gtype);
end  

if ok1*ok2 > 0
  if s > 0
    if ptype == 0
      cl = error_message('Plotting the fields. A pcolor plot.');
    elseif ptype == 1
      cl = error_message('Plotting the fields. An arrow plot.');
    end
  end

  Window = figure(4);
  clf;
  set(Window,'name','MULTIPAC - Field Window'); 

  if gtype == 1
    load fields
    load geodata.n
    gridcons = geodata(1,1);
    plot_cavity_fields(ptype,gridcons,s);
  elseif gtype == 2
    load geodata.n
    gridcons = geodata(1,1);
    plot_coupler_fields(ptype,gridcons,s);
    if s > 0 & ptype == 0
      cl = error_message(['Fields are plotted, top electric-walls, bottom ' ...
			  'magnetic-walls.']);
    end
  else
    load geodatal.n
    gridcons = geodatal(1,1);
    plot_window_fields(ptype,gridcons,s);
    if s > 0 & ptype == 0
      cl = error_message(['Fields are plotted, top electric-walls, bottom ' ...
			  'magnetic-walls.']);
    end
  end

  if s > 0
    cl = error_message('                                                  ');
  end
end  

if cl == 1
  clear_window;
end  
% ----------------------------------------------------------------------
