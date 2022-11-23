% Program save_data.m
% -------------------------------------------------------------------------
% Saves the input and outdata into file mpgui_data.tar.gz.
%
% -------------------------------------------------------------------------
% CALLS TO : none
% 12/12/00 : Pasi Yla-Oijala, Rolf Nevalinna Institute
% -------------------------------------------------------------------------

error_message('                    ');
error_message('Saving the files ...');

load fieldparam
if fieldparam(1) == 1
  !tar cf mpgui_data.tar geodata.n param fieldparam counter_flevels.mat counter_initials.mat Acounter.mat Atcounter.mat Ccounter.mat Efcounter.mat fieldfile1.n model.mat mesh.mat kama_n.mat mpgui_inputs.mat

elseif fieldparam(1) == 2
  !tar cf mpgui_data.tar geodata.n param fieldparam counter_flevels.mat counter_initials.mat Acounter.mat Atcounter.mat Ccounter.mat Efcounter.mat fieldfileE.n fieldfileH.n model.mat mesh.mat kama_1.mat kama_2.mat alphas.mat mpgui_inputs.mat

elseif fieldparam(1) == 3
  !tar cf mpgui_data.tar geodatar.n geodataw.n geodatal.n param fieldparam counter_flevels.mat counter_initialsl.mat counter_initialsr.mat Acounterr.mat Acounterl.mat Atcounterr.mat Atcounterl.mat Ccounterr.mat Ccounterl.mat Efcounterr.mat Efcounterl.mat fieldfileE.n fieldfileH.n model.mat mesh.mat kama_1.mat kama_2.mat alphas.mat mpgui_inputs.mat
end

!gzip -9 mpgui_data.tar 

error_message('Input and output fiels saved into file mpgui_data.tar.gz.');
error_message('To load the saved files, choose Load in menu File.');
% -------------------------------------------------------------------------


