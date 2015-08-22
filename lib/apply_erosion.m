function [ eimg ] = apply_erosion( img, shape, val, varargin )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if ~ischar(shape), error('shape must be a char'); end

if nargin < 4
    se = strel(shape,val);
    eimg = imerode(img,se);
end

end

