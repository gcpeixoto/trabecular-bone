function Eapp = calcEapp(f,A,Li,Lf)
%
% calcEapp Computes the apparent Young's modulus.
%
% input: 
%   f: sum of reaction forces over the surface A
%   A: surface area
%   Li: initial length
%   Lf: final length
%
% output:
%   Eapp = apparent Young' modulus

s = (Li - Lf)/Li;   % strain

Eapp = f/(s*A);     % Young

end

