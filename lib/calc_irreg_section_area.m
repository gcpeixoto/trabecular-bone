function [ As ] = calc_irreg_section_area( flag, tol, node, elem )
%calc_irreg_section_area Computes the surface area for a 
% set of elements inside a given region bounded by 'tol'
%
%   flag: '+x'

% size checking
if  size(node,2) > 3
    node = node(:,1:3);
else
    error('node must be 3 columns');
end

if size(elem,2) == 3    
    d = 3;
elseif size(elem,2) == 4
    d = 4;
else
    error('elem must be 3 or 4 columns');
end

A = [];

% irregular section area in relation to xmax
if strcmp(flag,'+x') == 1
    n = node(:,1);
    nm = max( node(:,1) );
    ve = find( abs( n - nm ) <= tol ); % nodes whose zcoord is inside tol
        
    for e = 1:size(elem,1)
        count = 0;
        for col = 1:d                        
            for j = 1:length(ve)
                if elem(e,col) == ve(j) % finding which elements contain the nodes
                count = count + 1; % how many special nodes 
                end
            end
        end
        if count >= 2 % criterion: if 2 or more nodes, save the element
            A = [A; elem(e,:)];            
        end
    end
    if ~isempty(A)
        As = cumsum( patch_area(A,node) ); % compute areas of the saved elements
        As = As(end); % total area
        fprintf('Total irregular section area (+x) = %8.2f \n',As);
    else
        As = 0;
        warning('Area is zero');    
    end
    
end

% irregular section area in relation to xmin
if strcmp(flag,'-x') == 1
    n = node(:,1);
    nm = min( node(:,1) );
    ve = find( abs( n - nm ) <= tol );    
        
    for e = 1:size(elem,1)
        count = 0;
        for col = 1:d                        
            for j = 1:length(ve)
                if elem(e,col) == ve(j)
                count = count + 1;
                end
            end
        end
        if count >= 2
            A = [A; elem(e,:)];            
        end
    end
    if ~isempty(A)
        As = cumsum( patch_area(A,node) );
        As = As(end);
        fprintf('Total irregular section area (-x) = %8.2f \n',As);
    else
        As = 0;
        warning('Area is zero');
    end
    
end

% irregular section area in relation to ymax
if strcmp(flag,'+y') == 1
    n = node(:,2);
    nm = max( node(:,2) );
    ve = find( abs( n - nm ) <= tol ); % nodes whose zcoord is inside tol
        
    for e = 1:size(elem,1)
        count = 0;
        for col = 1:d                        
            for j = 1:length(ve)
                if elem(e,col) == ve(j) % finding which elements contain the nodes
                count = count + 1; % how many special nodes 
                end
            end
        end
        if count >= 2 % criterion: if 2 or more nodes, save the element
            A = [A; elem(e,:)];            
        end
    end
    if ~isempty(A)
        As = cumsum( patch_area(A,node) ); % compute areas of the saved elements
        As = As(end); % total area
        fprintf('Total irregular section area (+y) = %8.2f \n',As);
    else
        As = 0;
        warning('Area is zero');    
    end
    
end

% irregular section area in relation to ymin
if strcmp(flag,'-y') == 1
    n = node(:,2);
    nm = min( node(:,2) );
    ve = find( abs( n - nm ) <= tol );    
        
    for e = 1:size(elem,1)
        count = 0;
        for col = 1:d                        
            for j = 1:length(ve)
                if elem(e,col) == ve(j)
                count = count + 1;
                end
            end
        end
        if count >= 2
            A = [A; elem(e,:)];            
        end
    end
    if ~isempty(A)
        As = cumsum( patch_area(A,node) );
        As = As(end);
        fprintf('Total irregular section area (-y) = %8.2f \n',As);
    else
        As = 0;
        warning('Area is zero');
    end
    
end

% irregular section area in relation to zmax
if strcmp(flag,'+z') == 1
    n = node(:,3);
    nm = max( node(:,3) );
    ve = find( abs( n - nm ) <= tol ); % nodes whose zcoord is inside tol
        
    for e = 1:size(elem,1)
        count = 0;
        for col = 1:d                        
            for j = 1:length(ve)
                if elem(e,col) == ve(j) % finding which elements contain the nodes
                count = count + 1; % how many special nodes 
                end
            end
        end
        if count >= 2 % criterion: if 2 or more nodes, save the element
            A = [A; elem(e,:)];            
        end
    end
    if ~isempty(A)
        As = cumsum( patch_area(A,node) ); % compute areas of the saved elements
        As = As(end); % total area
        fprintf('Total irregular section area (+z) = %8.2f \n',As);
    else
        As = 0;
        warning('Area is zero');    
    end
    
end

% irregular section area in relation to zmin
if strcmp(flag,'-z') == 1
    n = node(:,3);
    nm = min( node(:,3) );
    ve = find( abs( n - nm ) <= tol );    
        
    for e = 1:size(elem,1)
        count = 0;
        for col = 1:d                        
            for j = 1:length(ve)
                if elem(e,col) == ve(j)
                count = count + 1;
                end
            end
        end
        if count >= 2
            A = [A; elem(e,:)];            
        end
    end
    if ~isempty(A)
        As = cumsum( patch_area(A,node) );
        As = As(end);
        fprintf('Total irregular section area (-z) = %8.2f \n',As);
    else
        As = 0;
        warning('Area is zero');
    end
    
end
        

end

