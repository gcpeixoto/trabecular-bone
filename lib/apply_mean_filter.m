function [ fimg ] = apply_mean_filter( img, value )
% apply filter to image

filt = value*ones(3,3); % defining pixel neighbourhood for filter; 
fimg = imfilter(img,filt);

end

