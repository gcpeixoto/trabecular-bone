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
sample = 'z265'; % sample (image sequence)

ls_dir = dir( fullfile( img_dir,fmt,sample,strcat('*.',fmt) ) );        

%% SAVING
% output MSH
svmsh = fullfile(save_dir,'/msh');
opsvmsh = true; % optional to save msh
opsvstl = false;

% output FEB
modelName=fullfile(feb_dir,'boneMultiregion');

% output DAT
out = fullfile(pwd,'../dat/rev');

%% PARAMETER SETTINGS

% volumetric mesh adjust
%{


* reratio (RER) r/d ranges     
    RER = 0.612  : regular tetrahedron    
    RER > 2.0    : bad shaped tetrahedron
    (===> But... it seems this is not the def. of 'cgalv2m'... 
    Setting 'rer' returns error while calling binary 'cgalmesh.mexmaci64'.
    Why?
    
%}

nimg = 5; % number of images to parse
method = 'cgalmesh';% meshing method;
isovalue = [];      % level-set value \phi; understand the formation law
rad = 3;           % Delaunay's sphere radius (triangle maximum size)
ang = 30;           % miminum angle of a surface triangle
rer = 0.612;       % maximum radius-edge ratio
dist = 0.5;         % maximum deviation from specified isosurfaces 

%% IMAGE PROCESSING

n = size(ls_dir,1); 
im = imread(ls_dir(1).name);
vimg = zeros( size(im,1), size(im,2), nimg ); % allocating image array

afilt = true;
sigma = 0.5;
kexp = 0.3;
szel = 3;
for i = 1:nimg
    
    if (nimg > n); error('Number of images to load exceeded!'); end;
        
    im = imread( ls_dir(i).name );    % read         
    im = rgb2gray(im);                % convert to gray
    
    if afilt                     
        % see filter types
        ftg = fspecial('gaussian',szel,sigma);         
        ft = kexp*ones(szel,szel);        
        imft = imfilter(im,ft);               % filtered image                                  
        imftg = imfilter(im,ftg);               % filtered image   
        
        %figure
        %subplot(1,3,1), subimage(im)
        %subplot(1,3,2), subimage(imft)
        %subplot(1,3,3), subimage(imftg)
        %imshowpair(im,imft,'montage')         % original,filtered
        %close all
    end
    
    level = multithresh(im);          % threshold    
    im = imquantize(im,level);        % determine region labels 
            
    vimg(:,:,i) = im;                 % image with 1,2       
    
    level2 = multithresh(imft);          
    imft = imquantize(imft,level2);      
    vimgft(:,:,i) = imft;                
    
end

vimg = uint8( vimg ); % convert from 'logical' to 'uint8'
vimgft = uint8( vimgft );

%% TETRAHEDRAL MESH GENERATION

% mesher and options 
opt.radbound = rad;
opt.angbound = ang;
opt.reratio = rer; 
opt.distbound = dist;

maxvol1 = 1;
maxvol2 = 1;
maxvol = strcat('1=',num2str(maxvol1),':2=',num2str(maxvol2)); % 'label1=size1:label2:size2'
dofix = 1;
ix = 1:size(vimg,1);
iy = 1:size(vimg,2);
iz = 1:size(vimg,3);

[node,elem,face] = vol2mesh(vimg,ix,iy,iz,opt,maxvol,dofix,method);
[node2,elem2,face2] = vol2mesh(vimgft,ix,iy,iz,opt,maxvol,dofix,method);

% find the surface triangles 
face11=volface(elem(elem(:,end)==1,1:4));  % bone
face22=volface(elem(elem(:,end)==2,1:4));  % marrow
facep=[face11;face22];

% find the interface shared by the two surfaces
[fc,count1,count2]=unique(sort(facep')','rows');  % find the duplicates
bins=hist(count2,1:size(facep,1));        % the duplicated triangles are on the interface
cc=bins(count2);


% 1.6. Mesh processing and reparing
% [node,elem] = meshcheckrepair(node_out,elem_out); % not good for binary

%% Plot
if strcmp(plt,'y')    
    
    plotmesh(node,face)
    %figure
    %plotmesh(node,face11)
    %figure
    %plotmesh(node,face22)
    %figure
    %plotmesh(node(:,1:3),facep(cc==2,1:3))  % plot the interface
    
    figure
    plotmesh(node2,face2)
    
end

% [no,fc]=binsurface(vimg(1:10,1:10,1:10)==1);
% plotmesh(no,fc);
% figure
% [no2,fc2]=binsurface(vimg(1:10,1:10,1:10)==2);
% plotmesh(no2,fc2);

%% MSH saving

%{
It's necessary to open the mesh in Gmsh and save it again in case of 
importing it into PreView.
%}
if opsvmsh
    savemsh(node,elem, fullfile(svmsh,'boneREV.msh'), {'trabeculae','marrow'});
    savemsh(node2,elem2, fullfile(svmsh,'boneREV2.msh'), {'trabeculae','marrow'});
end

if opsvstl    
    savestl(node(:,1:3),elem(:,1:3), fullfile(svmsh,'boneREV.stl'), 'REV');
    savestl(node2(:,1:3),elem2(:,1:4), fullfile(svmsh,'boneREV2.stl'), 'REV');
end