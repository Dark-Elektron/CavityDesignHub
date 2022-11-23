% Function program plot_triplots(side,s)
% -----------------------------------------------------------------------
% To plot triplots.
%
% -----------------------------------------------------------------------
% CALLS TO : plot_triplot.m, plot_triplot_coupler.m, plot_triplot_win.m
%            check_geotype.m
% 09/05/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
% 13/12/00 :                 - check_geotype.m added
% -----------------------------------------------------------------------

function plot_triplots(side,s)

if nargin < 2
  s = 1;
end

ok3 = check_geotype(side);
if ok3 == 0
  return;
end

cl = 0;
ok = check_fieldparam;

if ok > 0
  load fieldparam;
  gtype = fieldparam(1);

  if gtype == 1
    plot_triplot(s);
  elseif gtype == 2
    plot_triplot_coupler(1,s);
  elseif gtype == 3
    if nargin < 1
      load side
    end
    if side == 'c'
      cl=error_message(['To plot the triplot in window, choose warm side ' ...
			'or cold side']);
      return;
    end
    save side side
    
    plot_triplot_win(1,side,s);
  end
end

if cl == 1
  clear_window;
end  
% ---------------------------------------------------------------
