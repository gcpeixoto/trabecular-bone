clear all; close all; clc;

%% Image directory

format = 'tif'; % 'jpg';

if strcmp(format,'jpg') == 1
    
    base_dir = '/Users/gustavo/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/jpg269';
    addpath(base_dir); % adds dir to MATLAB path; otherwise 'imread' will return error 
    ls_dir = dir( fullfile(base_dir,'*.jpg'));

elseif strcmp(format,'tif') == 1

    base_dir = '/Users/gustavo/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/tiff269';
    addpath(base_dir); % adds dir to MATLAB path; otherwise 'imread' will return error 
    ls_dir = dir( fullfile(base_dir,'*.tif'));
end


%% Output
% output STL
output_dir_stl = '/Users/gustavo/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/stl/';

% output MSH
output_dir_msh = '/Users/gustavo/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/msh/';


%% Image processing
n = size(ls_dir,1);
im = imread(ls_dir(1).name);
vimg = zeros( size(im,1), size(im,2), n ); % both
vimg_in = zeros( size(im,1), size(im,2), n ); % marrow
vimg_out = zeros( size(im,1), size(im,2), n ); % bone

prep = 'n'; % pre-processing methods
sizelim = 5; % integer as the maximum pixel size of a isolated region
iter = 50;   % number of iterations for mesh smoothing operation
filterval = 1.5;

nimg = 240; % number of images to parse
for i = 1:nimg
    if (nimg > n); error('Number of images to load exceeded!'); end;
    % read and convert
    im = imread( ls_dir(i).name );          
        
    % thresholding
    im = uint8(im);    
    level = multithresh(im);
    im = imquantize(im,level);
        
    %cc = bwconncomp(im); % stores struct from connected components    
    %im = labelmatrix(cc); % labelling image regions (more efficient than 'bwimage'        
    
    vimg_out(:,:,i) = im;    
    vimg_in(:,:,i) = imcomplement( vimg_out(:,:,i) ); % complementary image    
        
    % 1.9: volumetric image processing
    switch prep
        case 'di' 
            % remove isolated regions: maybe not useful for trabeculae
            vimg_out(:,:,i) = deislands3d(vimg_out(:,:,i),sizelim);    
            vimg_in(:,:,i) = deislands3d(vimg_in(:,:,i),sizelim);                        
        case 'sm'
            % perform a memory-limited 3D image smoothing
            vimg_out(:,:,i) = smoothbinvol(vimg_out(:,:,i),iter);    
            vimg_in(:,:,i) = smoothbinvol(vimg_in(:,:,i),iter);                    
        case 'mf'
            % apply mean filter
            vimg_out(:,:,i) = apply_mean_filter(vimg_out(:,:,i),filterval);    
            vimg_in(:,:,i) = apply_mean_filter(vimg_in(:,:,i),filterval);                    
    end
end

%% Fill holes in volumetric image

vimg_out = uint8( vimg_out ); % convert from 'logical' to 'uint8'
vimg_in = uint8( vimg_in ); % convert from 'logical' to 'uint8'

% 1.9: volumetric image processing
maxgap = 10; % maximum gap size for image closing
vimg_out = fillholes3d(vimg_out,maxgap);

%% 1.1 streamlined mesh generation: v2m

isovalue = []; % level-set value \phi; understand the formation law
% maxsizetri = 2; % surface triangle char. length; unit: pixel, voxel??
maxvoltet = 1200; % tetrahedra volume; unit^3: pixel, voxel??
clear opt;
opt.radbound = 1000; % Delaunay's sphere radius
%opt.distbound = 0.5; % maximum deviation from the specified isosurfaces
%opt.autoregion = 1; % if 1, saves the interior points for each closed surface component
method = 'cgalmesh';
[node_out,elem_out,face_out] = v2m(vimg_out,isovalue,opt,maxvoltet,method);


%fname = '/Users/gustavo/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/csv/nodes.csv';
%z_nodes = find_nodes(node_out,0.1,'+z');
%csvwrite(fname,z_nodes');

% 1.6. Mesh processing and reparing
% [node,elem] = meshcheckrepair(node_out,elem_out); % not good for binary

%% Plot
plotmesh(node_out,face_out)
%figure
%plotmesh(node_in,face_in)


%savestl(node_out,elem_out, fullfile(output_dir_stl,'teste.stl'),'trabeculae');

% BUG? method savemsh outputs a mesh which is not being imported by PreView; 
% PreView message Failed importing Gmsh fle.
% Failed finding EndElements.
savemsh(node_out,elem_out, fullfile(output_dir_msh,'trabeculae-rev.msh'));
