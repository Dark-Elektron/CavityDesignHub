% Function program plot_quaplots(side,s)
% -----------------------------------------------------------------------
% To plot quaplots.
%
% -----------------------------------------------------------------------
% CALLS TO : plot_quaplot.m, plot_quaplot_coupler.m, plot_quaplot_win.m
%            check_geotype.m
% 25/09/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
% 13/12/00 :                 - check_geotype.m added
% -----------------------------------------------------------------------

function plot_quaplots(side,s)

if nargin < 2
  s = 1;
end

ok3 = check_geotype(side);
if ok3 == 0
  return;
end

ok = check_fieldparam;
if ok > 0
  load fieldparam;
  gtype = fieldparam(1);

  if gtype == 1
    plot_quaplot(s);
    
  elseif gtype == 2
    plot_quaplot_coupler(1,s);
    
  elseif gtype == 3
    if nargin < 1
      load side
    end
    if side == 'c'
      cl=error_message(['To plot the quaplot in window, choose warm side '...
		     'or cold side']);
      return;      
    end
    save side side

    plot_quaplot_win(1,side,s);
  end
end  
% ---------------------------------------------------------------
