% Function program run_mesh_generator.m
% ---------------------------------------------------------------------------
% Generates the mesh.
%
% ---------------------------------------------------------------------------
% CALLS TO : clear_window.m, test_geometry.m, check_geometry.m, 
%            error_message.m, make_model_cavity.m, make_model_coupler.m, 
%            make_model_window.m, plot_mesh.m
% 12/04/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% ---------------------------------------------------------------------------

function run_mesh_generator

clear_window;

okt = test_geometry(0);
if okt == 0
  return;
end

ok = check_geometry;

if ok == 1
  load fieldparam
  gtype   = fieldparam(1);
  if gtype == 1
    cl = error_message('Generating the mesh. A cavity.');
  elseif gtype == 2
    cl = error_message('Generating the mesh. A coupler.');
  elseif gtype == 3
    cl = error_message('Generating the mesh. A window.');
  end

  releps  = fieldparam(3);
  gridcons= fieldparam(4);

  if gtype == 1
    make_model_cavity(gridcons,releps);
  elseif gtype == 2
%    error_message('A coupler.');
    make_model_coupler(gridcons,releps);
  elseif gtype == 3
%    error_message('A window.');
    make_model_window(gridcons,releps);
  end  

  % plot the mesh  
  figure(3)
  plot_mesh;

  cl = error_message('Mesh generator finished.');
  cl = error_message('                        ');
end  
% --------------------------------------------------------------------------

