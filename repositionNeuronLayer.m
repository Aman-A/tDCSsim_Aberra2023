function repositionNeuronLayer(cell_id,cell_layer,layers_or_layer_set_num,...
                                            nrn_pop_name,nrn_model_ver,varargin)
%REPOSITIONNEURONLAYER Run intersection detection and repositioning algorithms
%on one cell model throughout layer of specified NeuronPop
%  
%   Inputs 
%   ------ 
%   cell_id : integer
%             number indicating which cell model to simulate
%             (cell_model_name = cellModelNames(cell_id))
%   layer_num : integer
%               number indicating index of layer within layer set in which
%               to simulate this cell
%   layers_or_layer_set_num : struct or integer
%                             either layer structure or layer_set_num, 
%                             used to load layer struture to populate with 
%                             neurons
%   nrn_pop_name : string
%                  name of neuron population (set of neuron models in
%                  layer and their azimuthal orientations)
%   nrn_model_ver : string
%                   name corresponding to set of neuron model parameters,
%                   specified in outputMorphParams
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
if nargin == 0
   cell_id = 6;
   cell_layer = 2;
   layers_or_layer_set_num = 10; 
   nrn_pop_name = 'nrn_pop1';
   nrn_model_ver = 'maxH';    
end
in.par_on = 1; % optional, set to 1 to use parfor loop to parallelize
in.cell_models_file = 'cell_models';
in.ROI_expand_factor = 1; % expand box around MeshROI 
in.cellROI_expand_factor = 0.5; % expand box around cell morphology to crop GM for intersection detection
in.use_end_points_morph = 1; % use branchpoints/terminations instead of full morphology
in.range_for_search = 0; %% default is to search by 1mm and else it is fixing by the layer boundary
in.mode = 'both'; % set the maximum searching distance to be 'both','cell_length','fixed', 'd2layer'
% if 'fixed', then the distance should be specified in 'dmax'
in.dmax = 1; % mm
in.repo_pos = 0; % 0 runs on all positions, otherwise input vector of position indices 
in.save = 1;
in.min_nearest_faces = 10;
in.shift_inside = 0; % Choose to move the outside neurons inside the original layers or not 
in.epsilon = 1e-5; % epsilon for extra distance to move the neurons inside the layer 
in = sl.in.processVarargin(in,varargin);
[mat_dir,data_dir] = addPaths_tDCSsim;
cell_model_name = cellModelNames(cell_id,'mat_dir',mat_dir,...
                        'cell_models_file',in.cell_models_file);                 
fprintf('Running intersection detection and repositioning algorithm on cell: %s\n',...
        cell_model_name); 
mode = in.mode;
dmax = in.dmax;
min_nearest_faces = in.min_nearest_faces;
%% Load base layers, NeuronPop, and MeshROI
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
% if (~strcmp(nrn_model_ver,'maxM'))
for i = [1,3,5] % expand ROI to include tets just outside
    ROI(i:i+1) = ROI(i:i+1)+in.ROI_expand_factor*diff(ROI(i:i+1))*[-1,1];
end
[GrayMatter.vertices, GrayMatter.faces] = clipMeshVertices(GrayMatter.vertices,...
                                                             GrayMatter.faces,...
                                                             ROI);
% end
% Crop GM mesh
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
cell_data.C = cell_data.C*1e-3;
cell_origins0 = layers(cell_layer).cell_origins; % initial cell origins
cell_normals = layers(cell_layer).cell_normals; % cell normals (defined by cell placement layer elements)
% extract struct fields used below to avoid broadcasting 
cellROI_expand_factor = in.cellROI_expand_factor; 
num_pos = layers(cell_layer).num_elem; 

if in.repo_pos == 0
    repo_pos = 1:num_pos;
else
    repo_pos = in.repo_pos;
end
% Initialize new structures
cell_origins = cell_origins0; % new cell origins - initialize to start
phis = phis0; % new azimuthal rotations (deg)
intersected_gm_before = zeros(num_pos,1); % 1 or 0 if original position intersected with GM before repositioning
intersected_gm_after = zeros(num_pos,1); % 1 or 0 if still intersects with GM after repositioning
shift_distances = zeros(num_pos,1); 
intersect_lower_boundary_before = zeros(num_pos,1);
intersect_lower_boundary_after = zeros(num_pos,1);
intersect_upper_boundary_before = zeros(num_pos,1);
intersect_upper_boundary_after = zeros(num_pos,1);
valid_rotations = zeros(num_pos,12); % boolean array for whether the rotation is valid or not
rotations_checked = zeros(num_pos,12); % the azimuthal rotation angle that was checked 
d2layer_before = zeros(num_pos,1);

if isfield(MeshROI,'full_surfacesV')
   gm = MeshROI.full_surfacesV.GrayMatter;
else
   gm = MeshROI.full_surfaces.GrayMatter;
end

inside_gm = inpolyhedron(gm,cell_origins0);

fprintf('Generated mesh %g cell bodies outside gray matter\n',sum(~inside_gm));
inside_gm_before = inside_gm;
numCPUs = feature('numCores'); % Get number of CPUs available
% Get lower boundary surface, either layer boundary or GM/WM boundary
layersP = loadLayers(layer_set_num,'opt','layersP');

if cell_layer < length(layers)
    lower_boundary_surf = layersP(2*cell_layer).surface;
else
    lower_boundary_surf = MeshROI.surfaces.WhiteMatter; % For placement of last layer (L6), the lower boundary is whitematter
end

if cell_layer == 1
    upper_boundary_surf = GrayMatter; % For placement of L1, the upper boundary is graymatter
else
    upper_boundary_surf = layersP(2*(cell_layer-1)).surface; % For other placement layer, the upper boundary is the lower boundary of previous layer
end



if numCPUs > 1 && in.par_on
    % PARALLEL LOOP
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
    
    if(isempty(dmax))
        dmax = 1;
        disp("optional argument 'dmax' should be specified if 'fixed' is passed in, 1mm is recommended");
    end
    
    
    parfor (i = repo_pos,poolobj.NumWorkers)
       % Get initial position
        cell_origin0i = cell_origins0(i,:); % initial cell origin of position i
        cell_normal0i = cell_normals(i,:); % initial cell normal
        phi0i = phis0(i); % initial azimuthal rotation (deg 0-360)
        
        Cposi = placeCell(cell_origin0i,cell_normal0i,C,phi0i); % get cell coordinates at this position
        % Run intersection detection with GM
        [found_intersectioni,gm_facesi,gm_vertsi] = ...
                                checkCellMeshIntersection(GrayMatter,Cposi,...                                                        
                                                        parent_inds,...
                                                        cellROI_expand_factor); 
                                                                          
        [~, ~,~,upper_distance] = checkInsideOriginalLayer(Cposi,cell_origin0i,...
                                                            cell_normal0i, ...
                                                            upper_boundary_surf,...
                                                            lower_boundary_surf);
        % check if soma outside but intersection not detected, expand bounding box
        if ~inside_gm(i) && ~found_intersectioni            
            if isempty(upper_distance)
                upper_distance = 0.5;
            end
            cell_origin0i = cell_origin0i - cell_normal0i * upper_distance;
            Cposi = placeCell(cell_origin0i,cell_normal0i,C,phi0i); % get cell coordinates at this updated position
            [found_intersectioni,gm_facesi,gm_vertsi] = ...
                                checkCellMeshIntersection(GrayMatter,Cposi,...                                                        
                                                          parent_inds,...
                                                          cellROI_expand_factor); 
        end
        
        % check lower layer boundary 
        [lower_intersect, upper_intersect,lower_distance,~] = ...
                checkInsideOriginalLayer(Cposi,cell_origin0i,cell_normal0i, ...
                                    upper_boundary_surf,lower_boundary_surf);        
       
        if isempty(lower_distance)
            lower_distance = 0;
        end
        intersect_lower_boundary_before(i) = lower_intersect;
        intersect_upper_boundary_before(i) = upper_intersect;
        d2layer_before(i) = lower_distance;
        
        if found_intersectioni
            intersected_gm_before(i) = 1;
            outer_meshi = struct('faces',gm_facesi,'vertices',gm_vertsi);            
            % calculate the range for fixing
            
            if(strcmp('both',mode) || strcmp('d2layer',mode))
                % calculate the distance to layer boundary
                % crop the boundary mesh 
                [layer_faces,layer_vertices] = getMeshNearCell(lower_boundary_surf,...
                                                    Cposi,cellROI_expand_factor,...
                                            'min_nearby_faces',min_nearest_faces);
                partial_lower_boundary_surf = struct('faces',layer_faces,...
                                                'vertices',layer_vertices);
                inverted_normal = -1*cell_normal0i;
                [flag,distance_to_boundary] = rayIntersectionMesh(cell_origin0i, ...
                            inverted_normal, partial_lower_boundary_surf);
                
                if (~flag)
                    % find the least distance of average of the three points
%                     vector_to_vertices = layer_vertices - cell_origin0i;
%                     distances = sqrt(sum(vector_to_vertices.^2,2));
%                     least_distances = mink(distances,3); % find the 3 least distances
%                     distance_to_boundary = mean(least_distances);
                    distance_to_boundary = 0;
                end
                
                if(strcmp('both',mode))  
                    range_for_search = min(distance_to_boundary,1); % compare 1mm with the distance to boundary
                else
                    range_for_search = distance_to_boundary;
                end
                
                [cell_origin_newi,phi_newi,fixed_intersectioni] = ...
                fixCellMeshIntersection_BinarySearch(cell_origin0i,cell_normal0i,...
                                                 C,phi0i,parent_inds,outer_meshi,...
                                                 'range',range_for_search);  
            elseif(strcmp('cell_length',mode))
                % just call BinarySearch function 
                % the default would just take the cell length for searching
                [cell_origin_newi,phi_newi,fixed_intersectioni] = ...
                    fixCellMeshIntersection_BinarySearch(cell_origin0i,cell_normal0i,...
                                                 C,phi0i,parent_inds,outer_meshi); 
            elseif(strcmp('fixed',mode))
                %if fixed, just use 1mm 
                [cell_origin_newi,phi_newi,fixed_intersectioni] = ...
                fixCellMeshIntersection_BinarySearch(cell_origin0i,cell_normal0i,...
                                     C,phi0i,parent_inds,outer_meshi,'range',dmax); 
            else
                error("unrecognized mode, possible option is 'fixed', 'both','cell_length','d2layer'")
            end
            
             cell_origins(i,:) = cell_origin_newi;
             phis(i) = phi_newi; 
             intersected_gm_after(i) = ~fixed_intersectioni;
             shift_distances(i) = vmag(cell_origin_newi - cell_origin0i);
             C_newposi = placeCell(cell_origin_newi,cell_normal0i,C,phi_newi); % get cell coordinates at this position
                          

             [lower_intersect, upper_intersect,lower_distance, upper_distance] = ...
                 checkInsideOriginalLayer(C_newposi,cell_origin_newi,...
                                         cell_normal0i,upper_boundary_surf,...
                                         lower_boundary_surf);
        
             intersect_lower_boundary_after(i) = lower_intersect;
             intersect_upper_boundary_after(i) = upper_intersect;

             % check for valid azimuthal rotations
             [valid_rotation,rotate_degrees] = checkValidRotations(cell_origins(i,:),cell_normals(i,:),...
                             C,phis(i),parent_inds,outer_meshi) 
             valid_rotations(i,:) = valid_rotation;
             rotations_checked(i,:) = rotate_degrees;
         
        else
            cell_origins(i,:) = cell_origins0(i,:);
            phis(i) = phis0(i);
            intersected_gm_after(i) = found_intersectioni;
            shift_distances(i) = 0; 
            intersect_lower_boundary_after(i) = lower_intersect;
            intersect_upper_boundary_after(i) = upper_intersect;
            
            outer_meshi = struct('faces',gm_facesi,'vertices',gm_vertsi);
            % check for valid azimuthal rotations
            [valid_rotation,rotate_degrees] = checkValidRotations(cell_origins(i,:),cell_normals(i,:),...
                             C,phis(i),parent_inds,outer_meshi) 
            valid_rotations(i,:) = valid_rotation;
            rotations_checked(i,:) = rotate_degrees;
            
        end 
                
        fprintf(['Finished position %g, found_intersection before = %g, after = %g ' ... 
            'Intersect lower boundary %g, intersect lower boundary after %g\n'],...
                i,intersected_gm_before(i),intersected_gm_after(i),intersect_lower_boundary_before(i),intersect_lower_boundary_after(i)); 
    end
    total_time = toc;
    delete(poolobj);
    if exist(pc_storage_dir,'dir')
        rmdir(pc_storage_dir,'s'); % delete parfor temporary files
    end
else
    % SERIAL LOOP - 1 cpu or nan (if env variable doesn't exist)
    tic
    
    for i = repo_pos
       % Get initial position
        cell_origin0i = cell_origins0(i,:); % initial cell origin of position i
        cell_normal0i = cell_normals(i,:); % initial cell normal
        phi0i = phis0(i); % initial azimuthal rotation (deg 0-360)
        
        Cposi = placeCell(cell_origin0i,cell_normal0i,C,phi0i); % get cell coordinates at this position
        % Run intersection detection with GM
        [found_intersectioni,gm_facesi,gm_vertsi] = ...
                            checkCellMeshIntersection(GrayMatter,Cposi,...                                                        
                                                      parent_inds,...
                                                      cellROI_expand_factor);         
        [~, ~,~,upper_distance] = checkInsideOriginalLayer(Cposi,cell_origin0i,...
                                                            cell_normal0i, ...
                                                            upper_boundary_surf,...
                                                            lower_boundary_surf);
        
        % soma outside but intersection not detected, expand bounding box
        if ~inside_gm(i) && ~found_intersectioni 
            if isempty(upper_distance)
                upper_distance = 0.5;
            end
            cell_origin0i = cell_origin0i - cell_normal0i * upper_distance;
            Cposi = placeCell(cell_origin0i,cell_normal0i,C,phi0i); % get cell coordinates at this updated position
            [found_intersectioni,gm_facesi,gm_vertsi] = ...
                                checkCellMeshIntersection(GrayMatter,Cposi,...                                                        
                                                          parent_inds,...
                                                          cellROI_expand_factor); 
        end
        
        % check lower layer boundary 
        [lower_intersect, upper_intersect,lower_distance,~] = ...
                checkInsideOriginalLayer(Cposi,cell_origin0i,cell_normal0i, ...
                                        upper_boundary_surf,lower_boundary_surf);
        
        if isempty(lower_distance)
            lower_distance = 0;
        end
        intersect_lower_boundary_before(i) = lower_intersect;
        intersect_upper_boundary_before(i) = upper_intersect;
        d2layer_before(i) = lower_distance;
        
        if found_intersectioni
            intersected_gm_before(i) = 1;
            outer_meshi = struct('faces',gm_facesi,'vertices',gm_vertsi);            
            % calculate the range for fixing
            
            if(strcmp('both',mode) || strcmp('d2layer',mode))
                % calculate the distance to layer boundary
                % crop the boundary mesh 
                [layer_faces,layer_vertices] = getMeshNearCell(lower_boundary_surf,...
                                                Cposi,cellROI_expand_factor,...
                                                'min_nearby_faces',min_nearest_faces);
                partial_lower_boundary_surf = struct('faces',layer_faces,...
                                                    'vertices',layer_vertices);
                inverted_normal = -1*cell_normal0i;
                [flag,distance_to_boundary] = rayIntersectionMesh(cell_origin0i, ...
                            inverted_normal, partial_lower_boundary_surf);
                
                if (~flag)
                    % find the least distance of average of the three points
%                     vector_to_vertices = layer_vertices - cell_origin0i;
%                     distances = sqrt(sum(vector_to_vertices.^2,2));
%                     least_distances = mink(distances,3); % find the 3 least distances
%                     distance_to_boundary = mean(least_distances);
                    distance_to_boundary = 0;
                end
                
                if(strcmp('both',mode))  
                    range_for_search = min(distance_to_boundary,1); % compare 1mm with the distance to boundary
                else
                    range_for_search = distance_to_boundary;
                end
                
                [cell_origin_newi,phi_newi,fixed_intersectioni] = ...
                fixCellMeshIntersection_BinarySearch(cell_origin0i,cell_normal0i,...
                                                 C,phi0i,parent_inds,outer_meshi,...
                                                 'range',range_for_search);  
            elseif(strcmp('cell_length',mode))
                % just call BinarySearch function 
                % the default would just take the cell length for searching
                [cell_origin_newi,phi_newi,fixed_intersectioni] = ...
                    fixCellMeshIntersection_BinarySearch(cell_origin0i,cell_normal0i,...
                                                 C,phi0i,parent_inds,outer_meshi); 
            elseif(strcmp('fixed',mode))
                %if fixed, just use 1mm 
                [cell_origin_newi,phi_newi,fixed_intersectioni] = ...
                fixCellMeshIntersection_BinarySearch(cell_origin0i,cell_normal0i,...
                                     C,phi0i,parent_inds,outer_meshi,'range',dmax); 
            else
                error("unrecognized mode, possible option is 'fixed', 'both','cell_length','d2layer'")
            end
            
             cell_origins(i,:) = cell_origin_newi;
             phis(i) = phi_newi; 
             intersected_gm_after(i) = ~fixed_intersectioni;
             shift_distances(i) = vmag(cell_origin_newi - cell_origin0i);
             C_newposi = placeCell(cell_origin_newi,cell_normal0i,C,phi_newi); % get cell coordinates at this position
             
             [lower_intersect, upper_intersect,lower_distance, upper_distance] = ...
                    checkInsideOriginalLayer(C_newposi,cell_origin_newi,...
                                            cell_normal0i,upper_boundary_surf,...
                                            lower_boundary_surf);
        
             intersect_lower_boundary_after(i) = lower_intersect;
             intersect_upper_boundary_after(i) = upper_intersect;
             
            % check the valid azimuthal rotations
            [valid_rotation,rotate_degrees] = checkValidRotations(cell_origins(i,:),...
                                                            cell_normals(i,:),...
                                                            C,phis(i),...
                                                            parent_inds,outer_meshi);
            valid_rotations(i,:) = valid_rotation;
            rotations_checked(i,:) = rotate_degrees;
             
        else
            cell_origins(i,:) = cell_origins0(i,:);
            phis(i) = phis0(i);
            intersected_gm_after(i) = found_intersectioni;
            shift_distances(i) = 0; 
            intersect_lower_boundary_after(i) = lower_intersect;
            intersect_upper_boundary_after(i) = upper_intersect;
            
            outer_meshi = struct('faces',gm_facesi,'vertices',gm_vertsi);
            
            % check for valid azimuthal rotations
            [valid_rotation,rotate_degrees] = checkValidRotations(cell_origins(i,:),cell_normals(i,:),...
                             C,phis(i),parent_inds,outer_meshi);
            valid_rotations(i,:) = valid_rotation;
            rotations_checked(i,:) = rotate_degrees;
        end   
        
        fprintf(['Finished position %g: found_intersection before = %g, after = %g, ' ... 
            'Intersect lower boundary before = %g, intersect lower boundary after %g\n'],...
                i,intersected_gm_before(i),intersected_gm_after(i),intersect_lower_boundary_before(i),intersect_lower_boundary_after(i)); 
    end
    total_time = toc; 
end
fprintf('Finished loop in %f sec\n',total_time); 
fprintf('Found %g intersections, %g remaining after fixing (%.2f %% fixed)\n',...
        sum(intersected_gm_before),sum(intersected_gm_after),...
        100*(1-sum(intersected_gm_after)/sum(intersected_gm_before)));
fprintf('Mean shift distance = %.3f mm (+/- %.3f mm std)\n',...
        mean(shift_distances(shift_distances>0)),std(shift_distances(shift_distances>0)));
fprintf('Range shift distances = [%.3f, %.3f] mm\n',min(shift_distances(shift_distances>0)),max(shift_distances)); 
% Save output   
save_data = struct(); 
save_data.settings = in; 
save_data.cell_id = cell_id; 
save_data.cell_layer = cell_layer; 
save_data.cell_origins0 = cell_origins0; % old origins
save_data.phis0 = phis0; % old azimuthal rotations
save_data.cell_origins = cell_origins; % new origins
save_data.cell_normals = cell_normals; % same normals
save_data.phis = phis; % new azimuthal rotations
save_data.intersected_gm_before = intersected_gm_before; % intersection with GM before fixing
save_data.intersected_gm_after = intersected_gm_after; % and after
save_data.shift_distances = shift_distances; 
save_data.total_time = total_time; % time to fix all cells in layer
save_data.intersect_lower_boundary_before = intersect_lower_boundary_before;
save_data.intersect_lower_boundary_after = intersect_lower_boundary_after;
save_data.d2layer_before = d2layer_before;
save_data.inside_gm_before = inside_gm_before;
save_data.intersect_upper_boundary_before = intersect_upper_boundary_before;
save_data.intersect_upper_boundary_after = intersect_upper_boundary_after;
save_data.valid_rotations  = valid_rotations;
save_data.rotations_checked = rotations_checked;

if isfield(MeshROI,'full_surfacesV')
   inside_gm_after = inpolyhedron(MeshROI.full_surfacesV.GrayMatter,cell_origins); % Check if all the cell origins are within the graymatter 
else
   inside_gm_after = inpolyhedron(MeshROI.full_surfaces.GrayMatter,cell_origins);
end
save_data.inside_gm_after = inside_gm_after; % check if all cells inside the gray matter mesh

fprintf('After fixing, %g cell bodies outside the gray matter \n',sum(~inside_gm_after)); 
% Path handling
nrn_pop_dir = fullfile(mat_dir,'output_data','layer_data',layers(1).mesh_name,...
                      [layers(1).mesh_name '_' layers(1).roi_name],...
                      layers(1).layer_set_name,nrn_model_ver); 
repos_dir = fullfile(nrn_pop_dir,[nrn_pop_name '_repos']); 
if ~exist(repos_dir,'dir')
   mkdir(repos_dir); 
   fprintf('Made folder for outputs: %s\n',repos_dir); 
end

if in.save == 1
    cell_file_name = sprintf('L%g_%s_repos_data_%s.mat', cell_layer, cell_model_name, mode); 
    cell_data_file = fullfile(repos_dir,cell_file_name);
    save(cell_data_file,'-STRUCT','save_data'); 
    fprintf('Saved data to %s\n',cell_data_file); 
end

end
    