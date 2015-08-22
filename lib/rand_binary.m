clear all; close all; clc;

nimg = 10;
npx = 30;

v = zeros(npx,npx,nimg);

%% arbitrary porous matrix 
for i = 1:nimg    
    A = randi(2,npx);
    A = A - 1;
    v(:,:,i) = A;
end

% conversion
v = uint8(v);

%% meshing 
isovalue = []; % level-set value \phi; understand the formation law
maxvoltet = 10; % tetrahedra volume; unit^3: pixel, voxel??
clear opt;
opt.radbound = 10; % Delaunay's sphere radius
method = 'cgalmesh';
[node,elem,face] = v2m(v,isovalue,opt,maxvoltet,method);

plotmesh(node,elem)