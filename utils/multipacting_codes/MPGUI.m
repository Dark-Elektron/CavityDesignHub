% =========================================================================
% ||                          MPGUI.m                                    ||
% =========================================================================
%
%           ================================================
%           |     MATLAB Graphical user interface for      |
%           | Multipacting simulation toolbox MultiPac 2.1 |
%           |       Rolf Nevanlinna Institute    2001      |
%           ================================================
%
%
% Graphical user interface (GUI) for the multipacting simulation toolbox 
% MultiPac 2.1 to be used with MATLAB 5.0 or later.
%
% ------------------------------------------------------------------------
% 14/03/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% 26/09/00 :                 - version 1.1
% 17/10/00 :                 - version 2.0
% 13/12/00 :                 - version 2.0 (with a reduced number of bugs)
% 14/07/00 :                 - version 2.1
% ------------------------------------------------------------------------

l = 'l'; r = 'r'; c = 'c';

CommandWindow = figure(1);
clf;
set(CommandWindow,'name','MULTIPAC - Command Window');
set(CommandWindow,'menubar','none');

% ------------- file menu --------------------
file = uimenu('Label','File');
uimenu(file,'Label','Restart','Callback','MPGUI');
uimenu(file,'Label','Clear Command Window','Callback','clear_window');
uimenu(file,'Label','Exit','Callback','close all');
uimenu(file,'Label','Save Input Window','Callback','save_parameters',...
       'Separator','on');
uimenu(file,'Label','Load Input Window','Callback','load_parameters');
uimenu(file,'Label','Save Data','Callback','save_data');
uimenu(file,'Label','Load Data','Callback','load_data');
uimenu(file,'Label','Delete Geometry','Callback','clear_geometry',...
       'Separator','on');
uimenu(file,'Label','Delete Fields','Callback','clear_fields');
uimenu(file,'Label','Delete Inputs','Callback','clear_inputs');
uimenu(file,'Label','Delete Outputs','Callback','clear_outputs');

% ------------- run menu --------------------
run = uimenu('Label','Run');
uimenu(run,'Label','MultiPac','Callback','run_multipac');

% EM fields
EMfields = uimenu(run,'Label','EM Fields','Separator','on');
uimenu(EMfields,'Label','Field Solver','Callback','run_field_solver');
uimenu(EMfields,'Label','Mesh Generator','Callback','run_mesh_generator',...
       'Separator','on');
uimenu(EMfields,'Label','Eigenvalues','Callback','eigen_solver');

% multipacting
multipacting = uimenu(run,'Label','Multipacting','Separator','on');
uimenu(multipacting,'Label','MP Analysis','Callback','run_mpanalysis');
counters = uimenu(multipacting,'Label','Counter Functions','Separator','on');
uimenu(counters,'Label','Cavity or Coupler','Callback',...
       'calculate_counts(''c'')');
uimenu(counters,'Label','Warm Side','Callback','calculate_counts(''l'')');
uimenu(counters,'Label','Cold Side','Callback','calculate_counts(''r'')');

% distance map
distancemap = uimenu(multipacting,'Label','Distance Map');
uimenu(distancemap,'Label','Cavity or Coupler','Callback',...
       'calculate_distance(1)');
uimenu(distancemap,'Label','Warm Side','Callback',...
       'calculate_distance_win(1,''l'')');
uimenu(distancemap,'Label','Cold Side','Callback',...
       'calculate_distance_win(1,''r'')');

% trajectory
trajectory = uimenu(multipacting,'Label','Trajectory');
uimenu(trajectory,'Label','Cavity or Coupler','Callback',...
       'calculate_trajectory(1,''c'')');
uimenu(trajectory,'Label','Warm Side','Callback',...
       'calculate_trajectory(1,''l'')');
uimenu(trajectory,'Label','Cold Side','Callback',...
       'calculate_trajectory(1,''r'')');

% demos
uimenu(run,'Label','Demo Cavity','Separator','on','Callback','demo_cavity');
uimenu(run,'Label','Demo Coupler','Callback','demo_coupler');
uimenu(run,'Label','Demo Window','Callback','demo_window');

fields = uimenu('Label','Fields');
uimenu(fields,'Label','Plot Mesh','Callback','plot_mesh');
uimenu(fields,'Label','Calculate Fields','Callback','calculate_fields',...
       'Separator','on');
plotfield = uimenu(fields,'Label','Plot FEM Fields','Separator','on');
uimenu(plotfield,'Label','Pcolor','Callback','plot_FEM_fields(0)');
uimenu(plotfield,'Label','Arrow','Callback','plot_FEM_fields(1)');

plotfieldmx = uimenu(fields,'Label','Plot Mixed Fields');
uimenu(plotfieldmx,'Label','Pcolor','Callback','plot_mixed_fields(0)');
uimenu(plotfieldmx,'Label','Arrow','Callback','plot_mixed_fields(1)');

inputs = uimenu('Label','Inputs');
uimenu(inputs,'Label','Create Inputs','Callback','create_inputs');
%plotinputs = uimenu(inputs,'Label','Plot Inputs');
uimenu(inputs,'Label','Plot Geometry','Callback','plot_geometry',...
       'Separator','on');
uimenu(inputs,'Label','Plot Initials','Callback','plot_initials(1)');
uimenu(inputs,'Label','Plot Secondary Yield','Callback',...
       'plot_secondary_yield');
uimenu(inputs,'Label','Plot Emission Angle','Callback','plot_initial_angle');

outputs = uimenu('Label','Outputs');
plcounter = uimenu(outputs,'Label','Plot Counter Function');
uimenu(plcounter,'Label','Cavity or Coupler','Callback','plot_counter');
uimenu(plcounter,'Label','Warm Side','Callback','plot_counter(''l'')');
uimenu(plcounter,'Label','Cold Side','Callback','plot_counter(''r'')');

pltotal = uimenu(outputs,'Label','Plot Total Counter');
uimenu(pltotal,'Label','Cavity or Coupler','Callback','plot_total');
uimenu(pltotal,'Label','Warm Side','Callback','plot_total(''l'')');
uimenu(pltotal,'Label','Cold Side','Callback','plot_total(''r'')');

plenhanced = uimenu(outputs,'Label','Plot Enhanced Counter');
uimenu(plenhanced,'Label','Cavity or Coupler','Callback','plot_enhanced');
uimenu(plenhanced,'Label','Warm Side','Callback','plot_enhanced(''l'')');
uimenu(plenhanced,'Label','Cold Side','Callback','plot_enhanced(''r'')');

plenergy = uimenu(outputs,'Label','Plot Impact Energy');
uimenu(plenergy,'Label','Cavity or Coupler','Callback','plot_energy');
uimenu(plenergy,'Label','Warm Side','Callback','plot_energy(''l'')');
uimenu(plenergy,'Label','Cold Side','Callback','plot_energy(''r'')');

triplot = uimenu(outputs,'Label','Plot TriPlot','Separator','on');
uimenu(triplot,'Label','Cavity or Coupler','Callback','plot_triplots(''c'')');
uimenu(triplot,'Label','Warm Side','Callback','plot_triplots(''l'')');
uimenu(triplot,'Label','Cold Side','Callback','plot_triplots(''r'')');

quaplot = uimenu(outputs,'Label','Plot QuaPlot');
uimenu(quaplot,'Label','Cavity or Coupler','Callback','plot_quaplots(''c'')');
uimenu(quaplot,'Label','Warm Side','Callback','plot_quaplots(''l'')');
uimenu(quaplot,'Label','Cold Side','Callback','plot_quaplots(''r'')');

distance = uimenu(outputs,'Label','Plot Distance Map','Separator','on');
uimenu(distance,'Label','Cavity or Coupler','Callback',...
       'plot_distance(0,''c'',1);');
uimenu(distance,'Label','Warm Side','Callback','plot_distance(0,''l'',1);');
uimenu(distance,'Label','Cold Side','Callback','plot_distance(0,''r'',1);');

trajectory = uimenu(outputs,'Label','Plot Trajectory');
uimenu(trajectory,'Label','Cavity or Coupler','Callback',...
       'plot_trajectory(''c'',1)');
uimenu(trajectory,'Label','Window','Callback','plot_trajectory(''l'',1)');

%uimenu('Label','Help','Callback',...
%       'clear_window;error_message(''Help is not implemented.'');');
%uimenu('Label','Help','Callback','online_help');
%HelpInputs = uimenu(Helps,'Label','Inputs','Callback','inputs_help');
%HelpOutputs= uimenu(Helps,'Label','Ouputs');
%HelpPlots  = uimenu(Helps,'Label','Plots');

% create a frame for the title
uicontrol(gcf,'Style','Frame','Units','Normalized','Position',...
	  [0.24 0.82 0.47 0.16],'BackgroundColor',[0.0 0.8 0.5]);
%	  [0.24 0.82 0.47 0.16],'BackgroundColor',[0.0 0.7 0.5]);
%	  [0.24 0.82 0.47 0.16],'BackgroundColor',[0.3 0.55 1.0]);
%	  [0.24 0.82 0.47 0.16],'BackgroundColor',[0.2 0.7 1.0]);

uicontrol('Style','Text','String','MATLAB Graphical user interface for',...
	  'Units','Normalized','Position',[0.28 0.91 0.42 0.05],...
	  'Max',2,'HorizontalAlignment','left','Fontsize',10,...
	  'BackgroundColor',[0.0 0.8 0.5]);
%	  'BackgroundColor',[0.0 0.7 0.5]);
%	  'BackgroundColor',[0.3 0.55 1.0]);
%	  'BackgroundColor',[0.2 0.7 1.0]);
uicontrol('Style','Text','String',...
	  'Multipacting simulation toolbox MultiPac 2.1',...
	  'Units','Normalized','Position',[0.25 0.87 0.45 0.05],...
	  'Max',2,'HorizontalAlignment','left','Fontsize',10,...
	  'BackgroundColor',[0.0 0.8 0.5]);
%	  'BackgroundColor',[0.0 0.7 0.5]);
%	  'BackgroundColor',[0.3 0.55 1.0]);
%	  'BackgroundColor',[0.2 0.7 1.0]);
uicontrol('Style','Text','String','      Rolf Nevanlinna Institute    2001',...
	  'Units','Normalized','Position',[0.25 0.83 0.45 0.05],...
	  'Max',2,'HorizontalAlignment','left','Fontsize',10,...
	  'BackgroundColor',[0.0 0.8 0.5]);
%	  'BackgroundColor',[0.0 0.7 0.5]);
%	  'BackgroundColor',[0.3 0.55 1.0]);
%	  'BackgroundColor',[0.2 0.7 1.0]);

% create a frame for the communication box
uicontrol(gcf,'Style','Frame','Units','Normalized',...
	  'Position',[0.09 0.09 0.82 0.70]);
uicontrol('Style','Text','String','','Units','Normalized','Position',...
	  [0.1 0.1 0.8 0.68],'Tag','CommandWin',...
          'Max',2,'HorizontalAlignment','left','Fontsize',10);
% -----------------------------------------------------------------------

%uimenu(run,'Label','Generate Initials','Callback','generate_initials(1)',...
%       'Separator','on');
