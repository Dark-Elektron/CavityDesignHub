% Function program U = pow2volt(P,R,Z)
% --------------------------------------------------------------------
% Transfors rf power to peak voltage for a coaxial coupler.
% INPUT  P : power(s) [W]
%        R : reflection coefficient
%        Z : impedance of the line
% --------------------------------------------------------------------
% CALLS TO : error_message.m
% 10/05/00 : Pasi Yla-Oijala - Rolf Nevanlinna Institute
% --------------------------------------------------------------------

function U = pow2volt(P,R,Z)

U = zeros(size(P));

if Z > 0
  U = (1+abs(R))*sqrt(2*Z*P);
else
  error_message('Impendace = 0. Voltage is not defined.');
end  
% --------------------------------------------------------------------

