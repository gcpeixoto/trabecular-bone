%% meshREVBone.m

%{
Metadata: 

Author: Gustavo Peixoto de Oliveira
Date: August, 2015
e-mail: gustavo.oliveira@ci.ufpb.br
Federal University of Paraiba

Description: mesher with Febio embedded structures intended to simulate
bone structures generated from image sequence reading. 

This script is intended to mesh both the trabecular (g) 
and bone marrow (p) region as follows:

REV (multiregion)
===

-------------------------------
|     |  |   /  /    /  \     |     g: grain region  (trabeculae) 1
|     \  \  /  /    /   |     |     p: porous region (marrow)     2
|      \  \/  /    /    |     |
|      |     /     |   /      |
|      /  g /      / g |   p  |
|     /     |     /    \      |
-------------------------------

----------------------------------------------------------------------
----------------------------------------------------------------------

Matlab toolbox and dependencies: 

* iso2mesh, by Qianqian Fang  @ <http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Home> 
* gibbon, by Kevin M. Moerman @ <http://www.gibboncode.org>

Add-ons 

* convertFeb2Xplt, by yskmt @Github

----------------------------------------------------------------------
----------------------------------------------------------------------

%}

%% DEFAULTS
clear all; close all; clc;

%% PATHS

% path
i2m_path = '/Users/gustavo/Programs/iso2mesh';
gib_path = 'Users/gustavo/Programs/gibbon';
img_dir = fullfile(pwd,'../img');
save_dir = fullfile(pwd,'../save');
feb_dir = fullfile(pwd,'../febio/rev');

addpath( genpath(i2m_path), genpath(gib_path), fullfile(pwd,'../lib'), genpath(img_dir));

% deletes old files
delete( fullfile(feb_dir,'/*.*') );
delete( fullfile(feb_dir,'/sim/*.*') );

%% DIARY
% creating and enabling diary
diary( fullfile(save_dir,'mesh_bone_multiregion.txt') );
diary on

disp('==== MESH REV BONE - MULTIREGION ====')

%% IMAGE DIRECTORY

fmt = 'jpg'; % image format 
sample = 'z269'; % sample (image sequence)

ls_dir = dir( fullfile( img_dir,fmt,sample,strcat('*.',fmt) ) );        

%% SAVING
% output MSH
svmsh = fullfile(save_dir,'/msh');
opsvmsh = true; % optional to save msh

% output FEB
modelName=fullfile(feb_dir,'boneMultiregion');

% output DAT
out = fullfile(pwd,'../dat/rev');

%% PARAMETER SETTINGS

iter = 1; % number of iterations for mesh smoothing operation
nimg = 240; % number of images to parse
maxgap = 3; % maximum gap size for image fill holes 
sizelim = 5; % integer as the maximum pixel size of a isolated region
filterval = 1.5;
prep = 'n'; % pre-processing methods

% volumetric mesh adjust
%{


* reratio (RER) r/d ranges     
    RER = 0.612  : regular tetrahedron    
    RER > 2.0    : bad shaped tetrahedron
    (===> But... it seems this is not the def. of 'cgalv2m'... 
    Setting 'rer' returns error while calling binary 'cgalmesh.mexmaci64'.
    Why?
    
%}
method = 'cgalmesh';% meshing method;
isovalue = [];      % level-set value \phi; understand the formation law
maxvoltet = 100;    % tetrahedra volume; unit^3: pixel, voxel??
rad = 100;          % Delaunay's sphere radius
ang = 30;           % miminum angle of a surface triangle
%rer = 0.7;          % maximum radius-edge ratio
dist = 1;           % maximum distance between the center of the surface 
                    % bounding circle and center of the element bounding sphere

% b.c. setting
top_tol = 0.5;
bot_tol = 0.5;
dpx = 0.0;
dpy = 0.0;
dpz = -1.0;

% how the boundary condition will be applied
%bctype = 'pressure';
bctype = 'displacement';


% plotting
viewAz = 65;
vieEl = 25;
plt = 'y'; % 'n' plot bone structure?

%% IMAGE PROCESSING

n = size(ls_dir,1); 
im = imread(ls_dir(1).name);
vimg = zeros( size(im,1), size(im,2), nimg ); % allocating image array

for i = 1:nimg
    if (nimg > n); error('Number of images to load exceeded!'); end;
    % read and convert
    im = imread( ls_dir(i).name );          
        
    % thresholding       
    im = rgb2gray(im); % convert to gray
    level = multithresh(im); % determining regions
    im = imquantize(im,level);
    
    vimg(:,:,i) = im;        
        
    % 1.9: volumetric image processing
    switch prep
        case 'di' 
            % remove isolated regions: maybe not useful for trabeculae
            vimg(:,:,i) = deislands3d(vimg(:,:,i),sizelim);               
        case 'sm'
            % perform a memory-limited 3D image smoothing
            vimg(:,:,i) = smoothbinvol(vimg(:,:,i),iter);                
        case 'mf'
            % apply mean filter
            vimg(:,:,i) = apply_mean_filter(vimg(:,:,i),filterval);    
    end
end

vimg = uint8( vimg ); % convert from 'logical' to 'uint8'

%% FILL HOLES
% 1.9: volumetric image processing

vimg = fillholes3d(vimg,maxgap);

%% TETRAHEDRAL MESH GENERATION

% mesher and options 
opt.radbound = rad;
opt.angbound = ang;
%opt.reratio = rer; % error, if enabled
%opt.distbound = dist;
[node,elem,face] = v2m(vimg,isovalue,opt,maxvoltet,method);

% 1.6. Mesh processing and reparing
% [node,elem] = meshcheckrepair(node_out,elem_out); % not good for binary

%% Plot
if strcmp(plt,'y')
    plotmesh(node,face)
end

%% MSH saving

%{
It's necessary to open the mesh in Gmsh and save it again in case of 
importing it into PreView.
%}
if opsvmsh
    savemsh(node,elem, fullfile(svmsh,'boneREV.msh'), {'trabeculae','marrow'});
end