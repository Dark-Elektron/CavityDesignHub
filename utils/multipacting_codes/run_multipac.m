% Function program run_multipac.m
% --------------------------------------------------------------------
% Carried out the entire MP analysis with the field, counter
% function, distance map and electron trajectory calculations.
%
% --------------------------------------------------------------------
% CALLS TO : clear_window.m, error_message.m, test_geometry.m, 
%            run_field_solver.m, run_mpanalysis.m
% 10/05/00 : Pasi Yla-Oijala, Rolf Nevanlinna Institute
% 17/10/00 :                - preprocessor removed
% ---------------------------------------------------------------------

function MP_analysis

clear_window;

error_message('---------- Running MultiPac 2.0 ----------');
error_message('                                          ');

pause(3);

ok = test_geometry(0);
if ok == 0
  return;
end

% solve the EM fields
run_field_solver(0);

% carry out the multipacting analysis
run_mpanalysis(0);

error_message('                                          ');
error_message('--------- MultiPac 2.0 completed ---------');
% ---------------------------------------------------------------------

