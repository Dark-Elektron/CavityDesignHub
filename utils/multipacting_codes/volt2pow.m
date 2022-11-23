% Function program P = volt2pow(U,R,Z)
% --------------------------------------------------------------------
% Transforms peak voltage to rf power.
% INPUT  U : voltage(s) [V]
%        R : reflection coefficient
%        Z : impedance of the line
% --------------------------------------------------------------------
% CALLS TO : error_message.m
% 10/05/00 : Pasi Ylä-Oijala - Rolf Nevanlinna Institute
% --------------------------------------------------------------------

function P = volt2pow(U,R,Z)

P = zeros(size(U));

if Z > 0
  P = (U.^2)/(2*Z*((1+abs(R))^2));
else
  error_message('Impendace = 0. Power is not defined.');
end  
% --------------------------------------------------------------------

