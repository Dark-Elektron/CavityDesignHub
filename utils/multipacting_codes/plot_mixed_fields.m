% Function program plot_mixed_fields(ptype,s)
% ----------------------------------------------------------------------
% Plots the mixed fields used in the multipacting analysis.
% INPUT  ptype : 0 - pcolor, 1 - arrow
% ----------------------------------------------------------------------
% CALLS TO : check_geodata.m, error_message.m, plot_fields.m, 
%            clear_window.m, check_fields.m
% 10/05/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
% 17/10/00 :                 - coefficient files removed
% 13/12/00 :                 - check_fields.m added
% ----------------------------------------------------------------------

function plot_mixed_fields(ptype,s)

if nargin < 3
  s = 1;
end

cl  = 0;
ok1 = check_geodata;
ok2 = check_fieldparam;
if ok2 > 0
  load fieldparam
  gtype = fieldparam(1);
  ok3 = check_fields(gtype);
else
  ok3 = 0;
end  

if ok1*ok2*ok3 > 0
  load fieldparam
  gtype = fieldparam(1);

  if gtype == 1
    cl = error_message('Mixed fields are not defined for a cavity.');      
  else
    if s > 0
      if ptype == 0
	cl = error_message('Plotting the mixed fields. A pcolor plot.');
      elseif ptype == 1
	cl = error_message('Plotting the mixed fields. An arrow plot.');
      end
    end

    mixed_waves(s);

    load mixfields
    nr  = length(r);
    nz  = length(z);
    ERr = zeros(nr,nz);
    EZr = ERr; BPr = ERr;
    ERi = zeros(nr,nz);
    EZi = ERi; BPi = ERi;

    ERr(:) = real(Er);
    EZr(:) = real(Ez);
    BPr(:) = real(Bp);

    ERi(:) = imag(Er);
    EZi(:) = imag(Ez);
    BPi(:) = imag(Bp);

    if gtype == 2
      load geodata.n
      gridcons = geodata(1,1);
    else
      load geodatal.n
      gridcons = geodatal(1,1);
    end

    Window = figure(4);
    clf;
    set(Window,'name','MULTIPAC - Field Window'); 

    plot_fields(ptype,gridcons,ERr,EZr,BPr,ERi,EZi,BPi,r,z);
      
    if s > 0 
      if ptype == 0
	cl = error_message(['Fields are plotted, top electric field, ' ...
			    'bottom magnetic field.']);
      end
      cl = error_message('                                                ');
    end
  end
end  

if cl == 1
  clear_window;
end
% ----------------------------------------------------------------------

