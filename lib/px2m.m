function y = px2m( x, n, varargin )
% px2m Convert from px to meter
%
%   If n = 1, x must be a value in px   (line)
%   If n = 2, x must be a value in px^2 (area)
%   If n = 3, x must be a value in px^3 (volume)
%
%   input:
%       x: number in pixels
%       n: dimension
%
%   output
%       y: number in meters

% assumed conversion: 1 px = 34 micrometers; Roque (2013)
m = (34.0*1e-6);

if nargin == 1 
    y = x*m;
else
    switch n 
        case 1
            y = x*m;    % value in m   (line)
        case 2
            y = x*m^2;  % value in m^2 (area)
        case 3
            y = x*m^3;  % value in m^3 (volume) 
    end
end

end

