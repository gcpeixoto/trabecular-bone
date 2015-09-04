function y = px2m( x, n, varargin)
% px2m Convert from px to meter
%
%   input:
%       y: number in pixels
%       n: dimension
%
%   output
%       x: number in meters

% 1 px = 34 micrometers; Roque (2013)
m = (34*1e-6);

if nargin == 1
    y = x*m;
else
    switch n 
        case 1
            y = x*m;      % value in m
        case 2
            y = x*m^2;  % value in m^2 (area)
        case 3
            y = x*m^3;  % value in m^3 (volume) 
    end
end

end

