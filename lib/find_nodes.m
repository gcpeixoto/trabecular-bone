function [ varargout ] = find_nodes( v, eps, coord )
% find nodes in the neighborhood of extrema points
% [vnodes, ~]: return array of node indexing
% [vnodes, max]: return array of node indexing + extremum coordinate
%

% checking
if size(v,2) < 3, error('v must be at least 3-columns'); end

% search extrema points
if strcmp(coord,'+x') == 1
    a = v(:,1);
    am = max(a);    
    vnodes = find( abs( a(:) - am ) <= eps );        

elseif strcmp(coord,'-x') == 1
    a = v(:,1);
    am = min(a);
    vnodes = find( abs( a(:) - am ) <= eps );        

elseif strcmp(coord,'+y') == 1
    a = v(:,2);
    am = max(a);
    vnodes = find( abs( a(:) - am ) <= eps );            

elseif strcmp(coord,'-y') == 1
    a = v(:,2);
    am = min(a);
    vnodes = find( abs( a(:) - am ) <= eps );        

elseif strcmp(coord,'+z') == 1
    a = v(:,3);
    am = max(a);    
    vnodes = find( abs( a(:) - am ) <= eps );        

elseif strcmp(coord,'-z') == 1
    a = v(:,3);
    am = min(a);    
    vnodes = find( abs( a(:) - am ) <= eps );        
end

% print output 
if isempty(vnodes)
    warning('array of nodes is empty.\n'); 
else
    len = length(vnodes);
    st = fprintf('nodes found: %d ',len);
    fprintf( strcat( st, fprintf( 'extremum = %f \n', am ) ) );
    i = 1;
    while (i <= len)        
        %fprintf( 'node: %d \n', vnodes(i) );
        i = i+1;
    end
end
    
switch nargout
    
    % only array of nodes
    case 1
    varargout{1} = vnodes;
    
    % array of nodes and extremum coordinate
    case 2
    varargout{1} = vnodes;
    varargout{2} = am;
end

end

