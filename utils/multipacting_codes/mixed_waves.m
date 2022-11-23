% Function program mixed_waves(s)
% ---------------------------------------------------------------------
% Generates mixed waves.
%
% ---------------------------------------------------------------------
% CALLS TO : clear_window.m, check_fieldparam.m, error_message.m, mix.m
% 08/05/00 : Pasi Ylä-Oijala - Rolf Nevalinna Institute
% 17/10/00 :                 - coefficient files removed
% ---------------------------------------------------------------------

function mixed_waves(s)

if nargin < 1
  s = 1;
end  

cl = 0;
ok1 = check_fieldparam;

if ok1 > 0
  load fieldparam
  gtype = fieldparam(1);
  if gtype == 1
    cl = error_message('Mixed waves are not defined for a cavity.');
    return;
  else
    R = fieldparam(5) + i*fieldparam(6);

    % mix the coefficient files
    [alp1,alp2] = mix(R);
    disp('Ref coeff         alpha1              alpha2')
    disp(['  ' num2str(R) '   ' num2str(alp1) '    ' num2str(alp2)])
    disp('                            ')

    load fields1
    Er1 = i*Er;
    Ez1 = i*Ez;
    Bp1 = 4e-7*pi*H;
    
    load fields2
    Er2 = i*Er;
    Ez2 = i*Ez;
    Bp2 = 4e-7*pi*H;
    r   = I2;
    z   = I1;
    
    Er = alp1*Er1 + alp2*Er2;
    Ez = alp1*Ez1 + alp2*Ez2;
    Bp = alp1*Bp1 + alp2*Bp2;
    
    save mixfields Er Ez Bp rr zz r z
    
    if gtype == 2
      if s > 1
	cl = error_message('                                     ');
	cl = error_message('Generating mixed waves for a coupler.');
      elseif s > 0
	cl = error_message('Generating mixed waves for a coupler.');
      end
      
    else
      if s > 0
	cl = error_message('Generating mixed waves for a window.');
      end
    end
    if s > 1
      cl = error_message('Generation of mixed waves finished.');
      cl = error_message('                                   ');
    end
  end
end

if cl == 1
  clear_window;
end  
% ----------------------------------------------------------
