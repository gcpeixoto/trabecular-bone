clear all; close all; clc;

%% PATHS
i2m_path = '/Users/gustavo/Programs/iso2mesh';
img_dir = fullfile(pwd,'../img');
save_dir = fullfile(pwd,'../save');
feb_dir = fullfile(pwd,'../febio');

addpath( genpath(i2m_path), fullfile(pwd,'../lib'), genpath(img_dir));

%% IMAGE DIRECTORY

fmt = 'jpg';

if strcmp(fmt,'jpg') == 1    
    bdir = fullfile(img_dir,'/jpg/z269');    
    ls_dir = dir( fullfile(bdir,'*.jpg'));
end

%% SAVING
% output MSH
svmsh = fullfile(save_dir,'/msh');

% output FEB
modelName=fullfile(feb_dir,'boneCompression');

%% PARAMETER SETTINGS

iter = 10; % number of iterations for mesh smoothing operation
nimg = 80; % number of images to parse
maxgap = 10; % maximum gap size for image fill holes 

% image smoothing
proc = 'smooth';

% volumetric mesh
isovalue = []; % level-set value \phi; understand the formation law
maxvoltet = 120; % tetrahedra volume; unit^3: pixel, voxel??
rdbound = 10; % Delaunay's sphere radius
method = 'cgalmesh'; % meshing method;

% b.c. setting
top_tol = 1e-1;
bot_tol = 1e-1;
dpx = 0.0;
dpy = 0.0;
dpz = -1.5;

% plotting
viewAz = 65;
vieEl = 25;
plt = 'y'; % 'n' plot bone structure?

%% IMAGE PROCESSING

n = size(ls_dir,1); 
im = imread(ls_dir(1).name);
vimg = zeros( size(im,1), size(im,2), n ); % image array
for i = 1:nimg
    if (nimg > n); error('Number of images to load exceeded!'); end;
    
    % read and convert
    im = imread( ls_dir(i).name );          
    
    level = graythresh(im); % uses Otsu's method to parse            
    im = im2bw(im,level); % convert to binary        
    vimg(:,:,i) = im;        
        
    if strcmp(proc,'smooth') == 1
        fprintf('Performing binary smoothing on image %d. \n',i);        
        % perform a memory-limited 3D image smoothing
        vimg(:,:,i) = smoothbinvol(vimg(:,:,i),iter);   
                
    end
    
end

vimg = uint8( vimg ); % convert from 'logical' to 'uint8'

%% FILL HOLES
% 1.9: volumetric image processing

vimg = fillholes3d(vimg,maxgap);

%% TETRAHEDRAL MESH GENERATION
% See Iso2mesh manual Section 1.1 streamlined mesh generation: v2m

clear opt;
opt.radbound = rdbound;
[node,elem,face] = v2m(vimg,isovalue,opt,maxvoltet,method);

% mesh reorientation
% Although I have to use this method, calling 'Invert'
% in PreView is still required. I don't know why...
elem(:,1:4) = meshreorient(node(:,1:3),elem(:,1:4));

%% FIND AND SET DISPLACEMENT BCs

% extrema nodes
[z_top,zmax] = find_nodes( node, top_tol, '+z');
[z_bot,zmin] = find_nodes( node, bot_tol, '-z');

% Define displacement magnitudes
if dpz > 0    
    warning('z-displacement B.C. magnitude was inverted from %f to %f ',dpz,-dpz);
    dpz = - dpz;    
end
displacementMagnitude=[dpx dpy dpz];

%% VERIFICATION PLOT
if strcmp(plt,'y'); 
    title('Trabecular bone structure');     
    view(viewAz,vieEl);
    camlight('headlight');
    lighting phong;    
    plotmesh(node,face)
end

%% MSH saving

% BUG? method savemsh outputs a mesh which is not being imported by PreView; 
% PreView message Failed importing Gmsh fle.
% Failed finding EndElements. 
%
% It's necessary to open the mesh in Gmsh and save it again in case of 
% importing it into PreView.

savemsh(node,elem, fullfile(svmsh,'trabecula.msh'), {'trabecula'});

%% CONSTRUCTING FEB MODEL (from Gibbon code example)

FEB_struct.febio_spec.version='2.0';
FEB_struct.Module.Type='solid';

% Defining file names
FEB_struct.run_filename=[modelName,'.feb']; %FEB file name
FEB_struct.run_logname=[modelName,'.txt']; %FEBio log file name

febMatID = ones( size( elem(:,1),1 ),1 );

%Creating FEB_struct
FEB_struct.Geometry.Nodes=node(:,1:3);
FEB_struct.Geometry.Elements={elem(:,1:4)}; %The element sets
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
FEB_struct.Control.Values={50,0.05,...
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
