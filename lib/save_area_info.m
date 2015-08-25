function [ ] = save_area_info( fname, SA, EA, TA, BA, AR, BR, GR )
% save_area_info Saves surface area information 
%
% input: parameters
% (SA) surface area
% (EA) effective area
% (+A) top area 
% (-A) bottom area 
% (+A)/(SA) alfa ratio 
% (-A)/(SA) beta ratio 
% (max)/(min) gamma ratio 
%
% output: write to file

fid = fopen( fname, 'w' );
fprintf( fid,' (SA) \t (EA) \t (+A) \t (-A) \t (+A)/(SA) \t (-A)/(SA) \t gamma_A \n');
fprintf( fid,' %8.2f \t %8.2f \t %8.2f \t %8.2f \t %8.2f \t %8.2f \t %8.2f \n',...
        SA, EA, TA, BA, AR, BR, GR );
fclose(fid);

end

