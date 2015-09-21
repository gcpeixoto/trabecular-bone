function [ nodebb ] = forcebbox( node, top_tol, bot_tol )
%   Forces z coordinates inside tol to be zmax or zmin 
%   to make a flat z at bounding box.

Z = node(:,3);

idmax = abs( max(Z) - Z ) <= top_tol;
idmin = abs( min(Z) + Z ) <= bot_tol;

Z(idmax) = max(Z);
Z(idmin) = min(Z);

nodebb = [ node(:,1) node(:,2) Z ];

end

