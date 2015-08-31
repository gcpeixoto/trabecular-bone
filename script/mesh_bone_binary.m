%% mesh_bone_binary.m

%{
Metadata: 

Author: Gustavo Peixoto de Oliveira
Date: August, 2015
e-mail: gustavo.oliveira@ci.ufpb.br
Federal University of Para?ba

Description: mesher with Febio embedded structures intended to simulate
bone structures generated from image sequence reading.

Matlab toolbox dependencies: 

* iso2mesh, by Qianqian Fang,
* gibbon, by Kevin M. Moerman

Troubleshooting notes:

* Sometimes, iso2mesh seems to create meshes with 'isolated nodes'. These
nodes may be observed in the node subplot created by calling 'plotmesh'. 
Generally, run the script and mesh again might correct this problem,
although I'm not sure if 'cgalmesh' uses some randomic process that 
doesn't produces the same mesh in each run. On the other hand, 
adjusting the following parameters may help:

- 'maxgap': increasing it produces more isolated nodes in the final mesh

- rdbound: helps to avoid isolated nodes;

* While running for material 'rigid body', check if the XML structure
in file .feb shows up 'mat=49' at 'id' field for the material. If so,
force by hand 'mat=1'. This is an error of the XML node builder 
that prevents PreView of reading the file. I don't know why this error
happens.

%}

%% DEFAULTS
clear all; close all; clc;

%% PATHS
i2m_path = '/Users/gustavo/Programs/iso2mesh';
gib_path = 'Users/gustavo/Programs/gibbon';
img_dir = fullfile(pwd,'../img');
save_dir = fullfile(pwd,'../save');
feb_dir = fullfile(pwd,'../febio');

addpath( genpath(i2m_path), genpath(gib_path), fullfile(pwd,'../lib'), genpath(img_dir));

% deletes old files
delete( fullfile(feb_dir,'/*.*') );
delete( fullfile(feb_dir,'/sim/*.*') );

%% DIARY
% creating and enabling diary
diary( fullfile(save_dir,'mesh_bone_binary_log.txt') );
diary on

disp('==== MESH_BONE_BINARY ====')

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

% output DAT
out = fullfile(pwd,'../dat/');


%% PARAMETER SETTINGS

iter = 1; % number of iterations for mesh smoothing operation
nimg = 4; % number of images to parse
maxgap = 1; % maximum gap size for image fill holes 

% image smoothing
proc = 'smooth';

% volumetric mesh
isovalue = []; % level-set value \phi; understand the formation law
maxvoltet = 0.5; % tetrahedra volume; unit^3: pixel, voxel??
rdbound = 2; % Delaunay's sphere radius
method = 'cgalmesh'; % meshing method;

% b.c. setting
top_tol = 0.3;
bot_tol = 0.22;
dpx = 0.0;
dpy = 0.0;
dpz = -0.1*nimg; % 'p% of height'; set 'scale factor' in Preview

% how the boundary condition will be applied
bctype = 'pressure';
%bctype = 'displacement';

% plotting
viewAz = 65;
vieEl = 25;
plt = 'y'; % 'n' plot bone structure?

%% FEBIO STRUCTURE CONTROL SETTINGS

%{ 

* When max_ups is set to 0, FEBio will use the Full-Newton method 
instead of the BFGS method. In other words, the stiffness matrix 
is reformed for every iteration. In this case it is recommended 
to increase the number of max_refs (to e.g. 50), since the default 
value might cause FEBio to terminate prematurely when convergence 
is slow.

* The lstol parameter controls the scaling of the vector direction 
obtained from the line search. A line search method is used to 
improve the convergence of the nonlinear (quasi-) Newton solution 
algorithm. After each BFGS or Newton iteration, this algorithm 
searches in the direction of the displacement increment for a solution 
that has less energy (that is, residual multiplied with the displacement 
increment) than the previous iteration.

See Section 3.5.1 of Febio Manual for default values

%}

ntime = 100; % number of iterations 
dt = 0.1; % initial time step
max_refs = 25; % max number of stiffness reformations (recalculation and factorization)
max_ups = 0; % max number of BFGS stiffness updates;
dtol = 0.001; % convergence tolerance on displacements
etol = 0.01; % convergence tolerance on energy
rtol = 0; % convergence tolerance on residual 
lstol = 0.9; % convergence tolerance on line search

% time step props
dtmin = 1e-4; % minimum time step size
dtmax = 0.05; % maximum time step size
max_re = 5; % max number of retries allowed per time step
opt_iter = 10; % optimal number of iterations
aggr = 1; % aggressiveness (?)

% material choice
mat = 'rigid body'; % not OK with constraint imposition
%mat = 'Mooney-Rivlin';
%mat = 'neo-Hookean'; % keep as default, although used for cartilage

if strcmp(mat,'Mooney-Rivlin');
    fprintf('----> Running for material: Mooney-Rivlin.\n');
    c1 = 1e-3; % coefficient of first invariant term
    c2 = 0; % coefficient of second invariant term
    k = c1*1e3; % bulk modulus

elseif strcmp(mat,'neo-Hookean');
    fprintf('----> Running for material: neo-Hookean.\n');
    E = 1e+9; % Young modulus
    nu = 0.3; % Poisson' ratio

elseif strcmp(mat,'rigid body');
    fprintf('----> Running for material: rigid body.\n');        
    % If the center_of_mass parameter is 
    % omitted, FEBio will calculate the center of mass automatically. 
    % In this case, a density must be specified. If the center_of_mass 
    % is defined the density parameter may be omitted. 
    % Omitting both will generate an error message.
    rho = 1030;     
    E = 1e+9; % Young modulus
    nu = 0.3; % Poisson' ratio    
else
    error('Model of material not specified. Stopping...');
end


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
        fprintf('----> Performing binary smoothing on image %d. \n',i);        
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

%% VERIFICATION PLOT
if strcmp(plt,'y'); 
    subplot(2,2,1)
    title('Trabecular bone structure');     
    view(viewAz,vieEl);
    camlight('headlight');
    lighting phong;    
    plotmesh(node,face)
    hold on  
    plotmesh(node,face,strcat('abs(z - max(z) ) <', num2str(top_tol)),'facecolor','y');
    plotmesh(node,face,strcat('abs(z - min(z) ) <', num2str(bot_tol)),'facecolor','b');
    
    subplot(2,2,2)
    title('Trabecular bone structure');     
    view(viewAz,vieEl);
    camlight('headlight');
    lighting phong;        
    plotmesh(node);    
    
    subplot(2,2,3)    
    title('Top BC faces');     
    view(viewAz,vieEl);
    camlight('headlight');
    lighting phong;    
    plotmesh(node,face,strcat('abs(z - max(z) ) <', num2str(top_tol)),'facecolor','y');
    
    subplot(2,2,4)    
    title('Bottom BC faces');     
    view(viewAz,vieEl);
    camlight('headlight');
    lighting phong;    
    plotmesh(node,face,strcat('abs(z - min(z) ) <', num2str(bot_tol)),'facecolor','b');
end

%% AREAS COMPUTATION
% computation of surface area through Gibbon
areas = patch_area( face(:,1:3),node(:,1:3) ); % area per element
As = cumsum(areas); 
As = As(end); % sum of all areas
fprintf('----> Total surface area = %8.2f \n',As);

% approximation for the irregular cross section surface areas
[top_area, elem_top] = calc_irreg_section_area('+z',top_tol,node,face(:,1:3)); 
[bot_area, elem_bot] = calc_irreg_section_area('-z',bot_tol,node,face(:,1:3)); 

effective_area = As - ( top_area + bot_area );
max_sec_area = max(top_area,bot_area);
min_sec_area = min(top_area,bot_area);

ratios_area = [ top_area/As bot_area/As max_sec_area/min_sec_area ];

%% FIND AND SET BCs

% extrema nodes
[z_top,zmax] = find_nodes( node, top_tol, '+z');
[z_bot,zmin] = find_nodes( node, bot_tol, '-z');

% Define displacement magnitude scale
% (?) It seems that the choice of direction (+ or -) has no effect
% on the simulated values. Depends on load curves only?
if dpz > 0    
    warning('z-displacement B.C. magnitude was inverted from %f to %f ',dpz,-dpz);
    dpz = - dpz;    
end
displacementMagnitude=[dpx dpy dpz];

% Define pressure magnitude
% (?) check if this is the real value wanted
pressureMagnitude = 1.0;

%% MSH saving

%{
BUG? method savemsh outputs a mesh which is not being imported by PreView; 
PreView message Failed importing Gmsh fle.
Failed finding EndElements. 

It's necessary to open the mesh in Gmsh and save it again in case of 
importing it into PreView.
%}
savemsh(node,elem, fullfile(svmsh,'trabecula.msh'), {'trabecula'});

%% CONSTRUCTING FEB MODEL (from Gibbon code example)

FEB_struct.febio_spec.version='2.0';
FEB_struct.Module.Type='solid';

% Defining file names
FEB_struct.run_filename=[modelName,'.feb']; %FEB file name
FEB_struct.run_logname=[modelName,'.txt']; %FEBio log file name

% Defining material ID
febMatID = ones( size( elem(:,1),1 ),1 );

%Creating FEB_struct
FEB_struct.Geometry.Nodes=node(:,1:3);
FEB_struct.Geometry.Elements={elem(:,1:4)}; %The element sets
FEB_struct.Geometry.ElementType={'tet4'}; %The element types
FEB_struct.Geometry.ElementMat={febMatID};
FEB_struct.Geometry.ElementsPartName={'Trabeculae'};

%Defining surfaces
FEB_struct.Geometry.Surface{1}.Set=elem_top;
FEB_struct.Geometry.Surface{1}.Type='tri3';
%FEB_struct.Geometry.Surface{1}.Name='Pressure_surface';


% DEFINING MATERIAL CONSTITUTIVE RELATION
if strcmp(mat,'Mooney-Rivlin')
    FEB_struct.Materials{1}.Type='Mooney-Rivlin';
    FEB_struct.Materials{1}.Name='bone';
    FEB_struct.Materials{1}.Properties={'c1','c2','k'};
    FEB_struct.Materials{1}.Values={c1,c2,k};

elseif strcmp(mat,'neo-Hookean')
    FEB_struct.Materials{1}.Type='neo-Hookean';
    FEB_struct.Materials{1}.Name='bone';
    FEB_struct.Materials{1}.Properties={'E','v'};
    FEB_struct.Materials{1}.Values={E,nu};

elseif strcmp(mat,'rigid body')
    FEB_struct.Materials{1}.Type='rigid body';
    FEB_struct.Materials{1}.Name='bone';
    FEB_struct.Materials{1}.Properties={'density','E','v'};
    FEB_struct.Materials{1}.Values={rho,E,nu};        
    
else
    error('Model of material not specified. Stopping...');
end
    
%Control section
FEB_struct.Control.AnalysisType='static';
FEB_struct.Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};
FEB_struct.Control.Values={ntime,dt,...
    max_refs,max_ups,...
    dtol,etol,rtol,lstol};
FEB_struct.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter','aggressiveness'};
FEB_struct.Control.TimeStepperValues={dtmin,dtmax,max_re,opt_iter,aggr};

%Defining node sets
FEB_struct.Geometry.NodeSet{1}.Set=z_bot;
FEB_struct.Geometry.NodeSet{1}.Name='fixedBC';
FEB_struct.Geometry.NodeSet{2}.Set=z_top;
FEB_struct.Geometry.NodeSet{2}.Name='prescribedBC';
   
if strcmp(bctype,'displacement')
    
    disp('----> Running bctype: displacement.')
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

    % adding BC information for rigid boundary
    if strcmp(mat,'rigid body') ~= 1 % if not rigid, prescribe 
        FEB_struct.Boundary.Fix{1}.bc='x';
        FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
        FEB_struct.Boundary.Fix{2}.bc='y';
        FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
        FEB_struct.Boundary.Fix{3}.bc='z';
        FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
    end

% Pressure loading
elseif strcmp(bctype,'pressure')   
    disp('----> Running bctype: pressure.')
    
    % zero displacement doesn't apply together with constraints
    % for rigid body, constraints are set below 
    if strcmp(mat,'rigid body') ~= 1 
        
        % zero displacement at bottom
        FEB_struct.Boundary.Fix{1}.bc='x';
        FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;

        FEB_struct.Boundary.Fix{2}.bc='y';
        FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;

        FEB_struct.Boundary.Fix{3}.bc='z';
        FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;        
        
    end

    % pressure condition at top
    FEB_struct.Loads.Surface_load{1}.Type='pressure';    
    FEB_struct.Loads.Surface_load{1}.lcPar='pressure';
    FEB_struct.Loads.Surface_load{1}.lcParValue=pressureMagnitude;
    FEB_struct.Loads.Surface_load{1}.lc=1;     
    FEB_struct.Loads.Surface_load{1}.Set=elem_top;    

end

if strcmp(mat,'rigid body') == 1
    
    % adds constraints to lower surface
    
    % DOF x
    FEB_struct.Constraints{1}.RigidId = '1';    
    FEB_struct.Constraints{1}.Fix{1}.bc='x';
    %FEB_struct.Constraints{1}.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
    FEB_struct.Constraints{1}.Fix{1}.Set=elem_bot;
    
    % DOF y
    FEB_struct.Constraints{1}.Fix{2}.bc='y';
    %FEB_struct.Constraints{1}.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
    FEB_struct.Constraints{1}.Fix{2}.Set=elem_bot;
    
    % DOF z
    FEB_struct.Constraints{1}.Fix{3}.bc='z';
    %FEB_struct.Constraints{1}.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
    FEB_struct.Constraints{1}.Fix{3}.Set=elem_bot;
    
    % DOF rotation x
    FEB_struct.Constraints{1}.Fix{4}.bc='Rx';
    %FEB_struct.Constraints{1}.Fix{4}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
    FEB_struct.Constraints{1}.Fix{4}.Set=elem_bot;
    
    % DOF rotation y
    FEB_struct.Constraints{1}.Fix{5}.bc='Ry';
    %FEB_struct.Constraints{1}.Fix{5}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
    FEB_struct.Constraints{1}.Fix{5}.Set=elem_bot;
    
    % DOF rotation z
    FEB_struct.Constraints{1}.Fix{6}.bc='Rz';
    %FEB_struct.Constraints{1}.Fix{6}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;    
    FEB_struct.Constraints{1}.Fix{6}.Set=elem_bot;
end

%Load curves 
% Fix inside PreView
% 3 load curves are required for each direction separately in displacement 
FEB_struct.LoadData.LoadCurves.id=1;
FEB_struct.LoadData.LoadCurves.type={'linear'};

if strcmp(bctype,'displacement')   
    FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 -1;]};
end
if strcmp(bctype,'pressure')   
    FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 1;]};
end


%Adding output requests
FEB_struct.Output.VarTypes={'displacement', 'relative volume','shell thickness',...
    'stress','reaction forces'};

%Specify log file output
run_node_output_name=[FEB_struct.run_filename(1:end-4),'_node_out.txt'];
FEB_struct.run_output_names={run_node_output_name};
FEB_struct.output_types={'node_data'};
FEB_struct.data_types={'ux;uy;uz'};

%% SAVING .FEB FILE

FEB_struct.disp_opt=0; %Display waitbars option
febStruct2febFile(FEB_struct);

%% SAVING MESH INFORMATION

save_area_info( fullfile(out,'areas.dat'), As, effective_area, ...
    top_area, bot_area, ratios_area(1), ratios_area(2), ratios_area(3) ); 

% close diary
diary off


