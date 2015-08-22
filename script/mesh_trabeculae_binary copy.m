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

% output CSV
output_dir_csv = '/Users/gustavo/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/csv/';

% output FEB
output_dir_feb = '/Users/gustavo/Dropbox/Cloud_Waldir_Gustavo/GPeixoto/feb/';
modelName=fullfile(output_dir_feb,'trabeculae-compression');


%% Image processing
n = size(ls_dir,1);
im = imread(ls_dir(1).name);
vimg_in = zeros( size(im,1), size(im,2), n ); % marrow
vimg_out = zeros( size(im,1), size(im,2), n ); % bone

iter_out = 10;   % number of iterations for mesh smoothing operation
iter_in = 1;

nimg = 80; % number of images to parse
for i = 1:nimg
    if (nimg > n); error('Number of images to load exceeded!'); end;
    % read and convert
    im = imread( ls_dir(i).name );          
    
    level = graythresh(im); % uses Otsu's method to parse            
    im = im2bw(im,level); % convert to binary        
    
    vimg_out(:,:,i) = im;    
    vimg_in(:,:,i) = imcomplement( vimg_out(:,:,i) ); % complementary image    
        
    % 1.9: volumetric image processing        
    % perform a memory-limited 3D image smoothing
    vimg_out(:,:,i) = smoothbinvol(vimg_out(:,:,i),iter_out);    
    %vimg_in(:,:,i) = smoothbinvol(vimg_in(:,:,i),iter_in);                    
    
end

%% Fill holes in volumetric image

vimg_out = uint8( vimg_out ); % convert from 'logical' to 'uint8'
vimg_in = uint8( vimg_in ); % convert from 'logical' to 'uint8'

% 1.9: volumetric image processing
maxgap = 10; % maximum gap size for image closing
vimg_out = fillholes3d(vimg_out,maxgap);

%% Tetrahedral mesh generation
% See Iso2mesh manual Section 1.1 streamlined mesh generation: v2m

isovalue = []; % level-set value \phi; understand the formation law
maxvoltet = 120; % tetrahedra volume; unit^3: pixel, voxel??
clear opt;
opt.radbound = 10; % Delaunay's sphere radius
method = 'cgalmesh';
[node_out,elem_out,face_out] = v2m(vimg_out,isovalue,opt,maxvoltet,method);

maxvoltet = 34; 
clear opt;
opt.radbound = 20;
[node_in,elem_in,face_in] = v2m(vimg_in,isovalue,opt,maxvoltet,method);

%% Mesh reorientation
% Although I have to use this method, calling 'Invert'
% in PreView is still required. I don't know why...

elem_out(:,1:4) = meshreorient(node_out(:,1:3),elem_out(:,1:4));

%% Find and set BC

% extrema nodes
[z_top,zmax] = find_nodes( node_out, 1e-2, '+z');
[z_bot,zmin] = find_nodes( node_out, 1e-2, '-z');

% Define displacement magnitudes
displacementMagnitude=[0 0 -1.5];


%% Verification Plot
subplot(1,2,1)
plotmesh(node_out,face_out)
%subplot(1,2,2)
%plotmesh(node_in,face_in)

%% MSH saving

% BUG? method savemsh outputs a mesh which is not being imported by PreView; 
% PreView message Failed importing Gmsh fle.
% Failed finding EndElements. 
%
% It's necessary to open the mesh in Gmsh and save it again in case of 
% importing it into PreView.

savemsh(node_out,elem_out, fullfile(output_dir_msh,'bone.msh'), {'bone'});
%savemsh(node_in,elem_in, fullfile(output_dir_msh,'marrow.msh'),{'marrow'});

%% CONSTRUCTING FEB MODEL (from Gibbon code example)

FEB_struct.febio_spec.version='2.0';
FEB_struct.Module.Type='solid';

% Defining file names
FEB_struct.run_filename=[modelName,'.feb']; %FEB file name
FEB_struct.run_logname=[modelName,'.txt']; %FEBio log file name

febMatID = ones(size(elem_out(:,1),1),1);

%Creating FEB_struct
FEB_struct.Geometry.Nodes=node_out(:,1:3);
FEB_struct.Geometry.Elements={elem_out(:,1:4)}; %The element sets
FEB_struct.Geometry.ElementType={'tet4'}; %The element types
FEB_struct.Geometry.ElementMat={febMatID};
FEB_struct.Geometry.ElementsPartName={'Trabeculae'};

% DEFINING MATERIALS
c1=1e-3;
k=c1*1e3;
FEB_struct.Materials{1}.Type='Mooney-Rivlin';
FEB_struct.Materials{1}.Name='cube_mat';
FEB_struct.Materials{1}.Properties={'c1','c2','k'};
FEB_struct.Materials{1}.Values={c1,0,k};

%Control section
FEB_struct.Control.AnalysisType='static';
FEB_struct.Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};
FEB_struct.Control.Values={100,0.05,...
    25,0,...
    0.001,0.01,0,0.9};
FEB_struct.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter','aggressiveness'};
FEB_struct.Control.TimeStepperValues={1e-4,0.05,5,10,1};

%Defining node sets
FEB_struct.Geometry.NodeSet{1}.Set=z_bot;
FEB_struct.Geometry.NodeSet{1}.Name='bcFixList';
FEB_struct.Geometry.NodeSet{2}.Set=z_top;
FEB_struct.Geometry.NodeSet{2}.Name='bcPrescribeList';

%Adding BC information
FEB_struct.Boundary.Fix{1}.bc='x';
FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
FEB_struct.Boundary.Fix{2}.bc='y';
FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
FEB_struct.Boundary.Fix{3}.bc='z';
FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;

FEB_struct.Boundary.Prescribe{1}.SetName=FEB_struct.Geometry.NodeSet{2}.Name;
FEB_struct.Boundary.Prescribe{1}.Scale=displacementMagnitude(1);
FEB_struct.Boundary.Prescribe{1}.bc='x';
FEB_struct.Boundary.Prescribe{1}.lc=1;
FEB_struct.Boundary.Prescribe{1}.Type='relative';
FEB_struct.Boundary.Prescribe{2}.SetName=FEB_struct.Geometry.NodeSet{2}.Name;
FEB_struct.Boundary.Prescribe{2}.Scale=displacementMagnitude(2);
FEB_struct.Boundary.Prescribe{2}.bc='y';
FEB_struct.Boundary.Prescribe{2}.lc=1;
FEB_struct.Boundary.Prescribe{2}.Type='relative';
FEB_struct.Boundary.Prescribe{3}.SetName=FEB_struct.Geometry.NodeSet{2}.Name;
FEB_struct.Boundary.Prescribe{3}.Scale=displacementMagnitude(3);
FEB_struct.Boundary.Prescribe{3}.bc='z';
FEB_struct.Boundary.Prescribe{3}.lc=1;
FEB_struct.Boundary.Prescribe{3}.Type='relative';

%Load curves
FEB_struct.LoadData.LoadCurves.id=1;
FEB_struct.LoadData.LoadCurves.type={'linear'};
FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 1;]};

%Adding output requests
FEB_struct.Output.VarTypes={'displacement','stress','relative volume','shell thickness'};

%Specify log file output
run_node_output_name=[FEB_struct.run_filename(1:end-4),'_node_out.txt'];
FEB_struct.run_output_names={run_node_output_name};
FEB_struct.output_types={'node_data'};
FEB_struct.data_types={'ux;uy;uz'};

%% SAVING .FEB FILE

FEB_struct.disp_opt=0; %Display waitbars option
febStruct2febFile(FEB_struct);
