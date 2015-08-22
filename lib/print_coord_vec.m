function [ coords ] = print_coord_vec( vnodes, vmesh )
% print coordinates of the set of nodes 'vnodes' contained in
% the global mesh nodes 'vmesh'

if ~isempty(vnodes) && ~isempty(vmesh)
    
    vnlen = length(vnodes);
    coords = zeros(vnlen,4);
    
    for i = 1:vnlen
        coords(i,1) = vnodes(i); % node ID
        coords(i,2:4) = vmesh( vnodes(i),1:3); % node coordinates (x,y,z)
    end
    format shortg;
    disp(coords);
        
end

end

