function CheckInsideLayer(cell_id,cell_layer,layers_or_layer_set_num,...
                                            nrn_pop_name,nrn_model_ver,varargin)
%%                                        
if nargin == 0
   cell_id = 1; 
   cell_layer = 1;
   layers_or_layer_set_num = 13; 
   nrn_pop_name = 'nrn_pop1';
   nrn_model_ver = 'maxH';    
end

if ismac || ispc
    par_on = 0;
else
    par_on = 1;  % set to 1 to use parfor loop
end
in.cell_models_file = 'cell_models';
in.ROI_expand_factor = 1; % expand box around MeshROI 
in.cellROI_expand_factor = 0.5; % expand box around cell morphology to crop GM for intersection detection
in.use_end_points_morph = 1; % use branchpoints/terminations instead of full morphology
in.range_for_search = 0; %% default is to search by 1mm and else it is fixing by the layer boundary
in.mode = 'check'; % set the maximum searching distance to be 'both','cell_length','fixed', 'd2layer'
% if 'fixed', then the distance should be specified in 'dmax'
in.dmax = 1;
in.plot_on = 0;
in.min_nearest_faces = 10;
in = sl.in.processVarargin(in,varargin);
[mat_dir,data_dir] = addPaths_tDCSsim;
cell_model_name = cellModelNames(cell_id,'mat_dir',mat_dir,...
                        'cell_models_file',in.cell_models_file);                 
fprintf('Running intersection detection and repositioning algorithm on cell: %s\n',...
        cell_model_name);
plot_on = in.plot_on;
mode = in.mode;
dmax = in.dmax;
min_nearest_faces = in.min_nearest_faces;
if isstruct(layers_or_layer_set_num)
   layers = layers_or_layer_set_num;
%    layer_set_num = layers(1).layer_set_num;
elseif isnumeric(layers_or_layer_set_num)
   layer_set_num = layers_or_layer_set_num;
   layers = loadLayers(layer_set_num);
end                    
NeuronPop0 = loadNeuronPop(layers,nrn_pop_name,nrn_model_ver);    
MeshROI = loadMeshROI(layers); 
% Get mesh structures
if isfield(MeshROI,'full_surfacesV')
    GrayMatter = MeshROI.full_surfacesV.GrayMatter; % get GM surface from FEM tetrahedral mesh if available
else
    GrayMatter = MeshROI.full_surfaces.GrayMatter;
end
ROI = MeshROI.ROI;
for i = [1,3,5] % expand ROI to include tets just outside
    ROI(i:i+1) = ROI(i:i+1)+in.ROI_expand_factor*diff(ROI(i:i+1))*[-1,1];
end
% Crop GM mesh
[GrayMatter.vertices, GrayMatter.faces] = clipMeshVertices(GrayMatter.vertices,...
                                                                 GrayMatter.faces,...
                                                                 ROI);
cell_ids = NeuronPop0.cell_ids;
% check that cell is placed in input layer and get azimuthal rotations
% (phis)
if any(cell_ids{cell_layer} == cell_id)
   % get phi from already generated NeuronPop
   phis = NeuronPop0.phis;
   phis0 = phis{cell_layer}{cell_ids{cell_layer} == cell_id};
else
    error('cell_id %g was not placed in layer %g in this NeuronPop (%s)',cell_id,cell_layer,nrn_pop_name);   
end
% Load cell data
cell_data = loadCellData(cell_id,nrn_model_ver); 
if in.use_end_points_morph
    % extract end points of cell morphology
   cell_data = extractCellSubMorph(cell_data,'all','sectype','end_points_only');
end
C = cell_data.C*1e-3; % cell coordinates in local coordinate system (convert to mm)
parent_inds = cell_data.parent_inds; % connectivity of cell coordinates
cell_origins0 = layers(cell_layer).cell_origins; % initial cell origins
cell_normals = layers(cell_layer).cell_normals; % cell normals (defined by cell placement layer elements)
% extract struct fields used below to avoid broadcasting 
cellROI_expand_factor = in.cellROI_expand_factor; 
num_pos = layers(cell_layer).num_elem; 

numCPUs = feature('numCores'); % Get number of CPUs available
% Get lower boundary surface, either layer boundary or GM/WM boundary
layersP = loadLayers(layer_set_num,'opt','layersP');
if cell_layer < length(layers)
    lower_boundary_surf = layersP(2*cell_layer).surface;
else
    lower_boundary_surf = MeshROI.surfaces.WhiteMatter; % For placement of last layer (L6), the lower boundary is whitematter
end
intersect_lower_boundary = zeros(num_pos,1); % if the cell normal intersects with the lower boundary 
intersect_lr_boundary_non = zeros(num_pos,1); % if the cell normal intersects with the lower boundary 
not_in_layer = zeros(num_pos,1);
tic; 
if numCPUs > 1 && par_on
     fprintf('Number of CPUs requested = %g\n',numCPUs);
    pc_storage_dir = fullfile(data_dir,'pc_storage',getenv('SLURM_JOB_ID'));
    if numCPUs > 4 % assume on cluster        
        mkdir(pc_storage_dir);
        pc = parcluster('local');
        pc.JobStorageLocation =  pc_storage_dir;
    else
       pc = parcluster('local'); % use default
    end        
    poolobj = parpool(pc,numCPUs);
    tic; 
    
    parfor (i = 1:num_pos,poolobj.NumWorkers)
            % Get initial position
            cell_origin0i = cell_origins0(i,:); % initial cell origin of position i
            cell_normal0i = cell_normals(i,:); % initial cell normal
            phi0i = phis0(i); % initial azimuthal rotation (deg 0-360)

            Cposi = placeCell(cell_origin0i,cell_normal0i,C,phi0i); % get cell coordinates at this position

            % calculate the distance to layer boundary
            % crop the boundary mesh 
            [layer_faces,layer_vertices] = getMeshNearCell(lower_boundary_surf,Cposi,cellROI_expand_factor);
            [current_layer_faces,current_layer_vertices] = getMeshNearCell(layers(cell_layer).surface,Cposi,cellROI_expand_factor);
            partial_lower_boundary_surf = struct('faces',layer_faces,'vertices',layer_vertices);
            partial_current_layer = struct('faces',current_layer_faces,'vertices',current_layer_vertices);
            
            inverted_normal = -1*cell_normal0i;
            % Check intersection of lower boundary with the inverted normals 
            [lower_intersect,~] = rayIntersectionMesh(cell_origin0i, inverted_normal, partial_lower_boundary_surf);
            intersect_lower_boundary(i) = lower_intersect;
            % Check intersection of lower boundary with the cell normals
            [lower_intersect_non_inverted,~] = rayIntersectionMesh(cell_origin0i, cell_normal0i, partial_lower_boundary_surf);
            intersect_lr_boundary_non(i) = lower_intersect_non_inverted;
            if(lower_intersect_non_inverted)
                disp(i);
                not_in_layer(i) = 1;
            end
    end
else
    for i = 1:num_pos
         cell_origin0i = cell_origins0(i,:); % initial cell origin of position i
            cell_normal0i = cell_normals(i,:); % initial cell normal
            phi0i = phis0(i); % initial azimuthal rotation (deg 0-360)

            Cposi = placeCell(cell_origin0i,cell_normal0i,C,phi0i); % get cell coordinates at this position

            % calculate the distance to layer boundary
            % crop the boundary mesh 
            [layer_faces,layer_vertices] = getMeshNearCell(lower_boundary_surf,Cposi,cellROI_expand_factor);
            [current_layer_faces,current_layer_vertices] = getMeshNearCell(layers(cell_layer).surface,Cposi,cellROI_expand_factor);
            partial_lower_boundary_surf = struct('faces',layer_faces,'vertices',layer_vertices);
            partial_current_layer = struct('faces',current_layer_faces,'vertices',current_layer_vertices);
            
            inverted_normal = -1*cell_normal0i;
            % Check intersection of lower boundary with the inverted normals 
            [lower_intersect,~] = rayIntersectionMesh(cell_origin0i, inverted_normal, partial_lower_boundary_surf);
            intersect_lower_boundary(i) = lower_intersect;
            % Check intersection of lower boundary with the cell normals
            [lower_intersect_non_inverted,~] = rayIntersectionMesh(cell_origin0i, cell_normal0i, partial_lower_boundary_surf);
            intersect_lr_boundary_non(i) = lower_intersect_non_inverted;
            if(lower_intersect_non_inverted)
                
                if (plot_on)
                    disp(i);
                    plb = plotDefPatch(partial_lower_boundary_surf); hold on;
                    plb.FaceAlpha = 0.7; 
                    pcl = plotDefPatch(partial_current_layer); hold on;
                    pcl.FaceAlpha = 0.7; 
                    pcl.FaceColor = 'r';
                    quiver3(cell_origin0i(1),cell_origin0i(2),cell_origin0i(3),cell_normal0i(1),cell_normal0i(2),cell_normal0i(3));
                    ax = gca; 
                    cell_data = struct('parent_inds',parent_inds); 
                    cell_data.C = placeCell(cell_origin0i,cell_normal0i,C,phi0i);
                    plotCellLines('cell_data',cell_data,'lw',2','vals',-10); 
                    ax_lims = axis; 
                    colormap(flipud([0,0,1;0,1,0;0.2*ones(1,3);1,0,0])); caxis([-10 20]); axis off; 
                    legend('lower boundary surface','current layer','cell normals');
                    drawnow;
                    [lower_intersect_non_inverted,~] = rayIntersectionMesh(cell_origin0i, cell_normal0i, partial_lower_boundary_surf);
                    intersect_lr_boundary_non(i) = lower_intersect_non_inverted;
                    figure();
                    quiver3(cell_origin0i(1),cell_origin0i(2),cell_origin0i(3),cell_normal0i(1),cell_normal0i(2),cell_normal0i(3));
                    plb = plotDefPatch(partial_lower_boundary_surf); hold on;
                    plb.FaceAlpha = 0.7; 
                end
               not_in_layer(i) = 1; 
            end
    end
end
disp(sum(not_in_layer));
save_data = struct(); 
save_data.settings = in; 
save_data.cell_id = cell_id; 
save_data.cell_layer = cell_layer; 
save_data.cell_origins0 = cell_origins0; % old origins
save_data.phis0 = phis0; % old azimuthal rotations
save_data.intersect_lower_boundary = intersect_lower_boundary;
save_data.intersect_lr_boundary_non = intersect_lr_boundary_non;
% Path handling
nrn_pop_dir = fullfile(mat_dir,'output_data','layer_data',layers(1).mesh_name,...
                      [layers(1).mesh_name '_' layers(1).roi_name],...
                      layers(1).layer_set_name,nrn_model_ver); 
repos_dir = fullfile(nrn_pop_dir,[nrn_pop_name '_repos']); 
if ~exist(repos_dir,'dir')
   mkdir(repos_dir); 
   fprintf('Made folder for outputs: %s\n',repos_dir); 
end
cell_file_name = sprintf('L%g_%s_repos_data_%s.mat',cell_layer,cell_model_name, mode); 
cell_data_file = fullfile(repos_dir,cell_file_name);
save(cell_data_file,'-STRUCT','save_data'); 
fprintf('Saved data to %s\n',cell_data_file); 

        
end