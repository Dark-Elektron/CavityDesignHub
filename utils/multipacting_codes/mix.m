% Function program [alpha1,alpha2] = mix(R)
% -----------------------------------------------------------------------
% Calculates coefficients for generating mixed waves by mixing two 
% solutions with electric and magnetic walls.
% INPUT  R : reflection coefficient
%      alpha1, alpha2 : coefficient for mixing, alpha1 for the solution
%              with E-walls and alpha2 for a solution with H-walls
% -----------------------------------------------------------------------
% CALLS TO : none
% xx/yy/00 : Marko Ukkola - Rolf Nevanlinna Institute
% 08/05/00 : Pasi Ylä-Oijala
% -----------------------------------------------------------------------

function [alpha1,alpha2] = mix(R)

load param

c     = 2.99792458e8;
freq  = param(1);
delta = c/(freq*4);
k     = 2 * pi * freq/c;
  
alpha2 = -2*(R + 1) ./ ((1+abs(R))*2*sqrt(-1)*sin(k*delta));
alpha1 = (2./(1+abs(R)) - alpha2 * exp(-sqrt(-1)*k*delta));

% save values for further use
alphas = [real(alpha1),imag(alpha1),real(alpha2),imag(alpha2)];
save alphas alphas -v4
% -----------------------------------------------------------------------
