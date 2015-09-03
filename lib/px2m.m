function y = px2m( x, n)
% px2m Convert from px to meter
%
%   input:
%       y: number in meters
%       n: dimension
%
%   output
%       x: number in pixels

% 1 voxel = 34 micrometers; then 1 px = cubic root of 34
px = (34*1e-6)^(1/3);

switch n 
    case 1
        y = x*px;      % value in m
    case 2
        y = x*(px)^2;  % value in m^2 (area)
    case 3
        y = x*(px)^3;  % value in m^3 (volume) 
end

end

