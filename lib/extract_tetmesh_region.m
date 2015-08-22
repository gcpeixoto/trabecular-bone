function [node,elem,face] = extract_tetmesh_region(n,el,fc,reg)
% extract mesh by region

if ( size(n,2) ~= 4 || size(el,2) ~=5 || size(fc,2) ~= 4 )
    error('Dimensions of input arrays do not determine to a tetrahedral mesh.');
end

col1 = []; 
col2 = [];
col3 = [];
col4 = [];
i = 1;
while i <= length(n)
    if n(i,4) == reg
        col1 = [col1; n(i,1) ];        
        col2 = [col2; n(i,2) ];
        col3 = [col3; n(i,3) ];
        col4 = [col4; n(i,4) ];
    end
    i = i + 1;
end

node = [col1 col2 col3 col4];

col1 = []; 
col2 = [];
col3 = [];
col4 = [];
col5 = [];
i = 1;
while i <= length(n)
    if el(i,5) == reg
        col1 = [col1; el(i,1) ];        
        col2 = [col2; el(i,2) ];
        col3 = [col3; el(i,3) ];
        col4 = [col4; el(i,4) ];        
        col5 = [col5; el(i,5) ];        
    end
    i = i + 1;
end

elem = [col1 col2 col3 col4 col5];

col1 = []; 
col2 = [];
col3 = [];
col4 = [];

i = 1;
while i <= length(n)
    if fc(i,4) == reg
        col1 = [col1; fc(i,1) ];        
        col2 = [col2; fc(i,2) ];
        col3 = [col3; fc(i,3) ];
        col4 = [col4; fc(i,4) ];              
    end
    i = i + 1;
end

face = [col1 col2 col3 col4];

end

