function compileRepositionedPop(layer_set_num,nrn_pop_name,nrn_model_ver,...
                                varargin)
%COMPILEREPOSITIONEDPOP Compile output of repositionNeuronLayer for input 
% neuron population 
%  
%   Inputs 
%   ------ 
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
if nargin == 0
   layer_set_num = 16; 
   nrn_pop_name = 'nrn_pop1';
   nrn_model_ver = 'maxM'; 
end
in.over_write = 0; % set to 1 to over write original NeuronPop file
in.cell_models_file = 'cell_models'; % default
in.mode = 'both'; % set the maximum searching distance to be 'both','cell_length','fixed', 'd2layer'
in = sl.in.processVarargin(in,varargin); 
mode = in.mode;
mat_dir = addPaths_tDCSsim;
layers = loadLayers(layer_set_num); % load original layers file
NeuronPop0 = loadNeuronPop(layers,nrn_pop_name,nrn_model_ver); 
nrn_pop_dir = fullfile(mat_dir,'output_data','layer_data',layers(1).mesh_name,...
                      [layers(1).mesh_name '_' layers(1).roi_name],...
                      layers(1).layer_set_name,nrn_model_ver); 
repos_dir = fullfile(nrn_pop_dir,[nrn_pop_name '_repos']); 
fprintf('shift mode: %s\n',mode); 
%% Loop through all cells and load new cell origins, azimuthal rotations, 
%  and stats
NeuronPop = NeuronPop0; % new NeuronPop
cell_ids = NeuronPop0.cell_ids; 
cell_ids_lin = [cell_ids{:}]; 
num_layers = length(cell_ids); 
% repositioned cell data
phis = cell(num_layers,1); 
cell_origins = cell(num_layers,1); 
cell_normals = cell(num_layers,1); 
% repositioning stats
intersected_gm_before = cell(num_layers,1); 
intersected_gm_after = cell(num_layers,1); 
stay_in_layer = cell(num_layers,1);
shift_distances = cell(num_layers,1); 
intersect_lower_boundary_before = cell(num_layers,1); 
intersect_lower_boundary_after = cell(num_layers,1); 
d2layer_before = cell(num_layers,1);
intersect_upper_boundary_before = cell(num_layers,1);
intersect_upper_boundary_after = cell(num_layers,1);
inside_gm_before = cell(num_layers,1);
inside_gm_after = cell(num_layers,1);

repositioned_cells = zeros(size(cell_ids_lin)); % binary, 1 if found repositioning data



for i = 1:num_layers   
    num_cells_in_layer = length(cell_ids{i}); 
    num_elem = layers(i).num_elem; 
    phis{i} = cell(1,num_cells_in_layer); 
    cell_origins{i} = cell(1,num_cells_in_layer); 
    cell_normals{i} = cell(1,num_cells_in_layer); 
    intersected_gm_before{i} = zeros(num_elem,num_cells_in_layer); 
    intersected_gm_after{i} = zeros(num_elem,num_cells_in_layer); 
    shift_distances{i} = zeros(num_elem,num_cells_in_layer); 
    stay_in_layer{i} = zeros(num_elem,num_cells_in_layer);
    intersect_lower_boundary_before{i} = zeros(num_elem,num_cells_in_layer);
    intersect_lower_boundary_after{i} = zeros(num_elem,num_cells_in_layer);
    d2layer_before{i} = zeros(num_elem,num_cells_in_layer);
    intersect_upper_boundary_before{i} = zeros(num_elem,num_cells_in_layer);
    intersect_upper_boundary_after{i} = zeros(num_elem,num_cells_in_layer);
    inside_gm_before{i} = zeros(num_elem,num_cells_in_layer);
    inside_gm_after{i} = zeros(num_elem,num_cells_in_layer);
    
   for j = 1:num_cells_in_layer
        cell_ij = cell_ids{i}(j);         
        cell_model_names_ij = cellModelNames(cell_ij,'mat_dir',mat_dir,...
                                            'cell_models_file',in.cell_models_file); 
        repos_file_name = sprintf('L%g_%s_repos_data_%s.mat',i,cell_model_names_ij,mode); 
        repos_file = fullfile(repos_dir,repos_file_name); 
        if exist(repos_file,'file')
            repos_data = load(repos_file); 
            cell_origins{i}{j} = repos_data.cell_origins; 
            cell_normals{i}{j} = repos_data.cell_normals; 
            phis{i}{j} = repos_data.phis; 
            intersected_gm_before{i}(:,j) = repos_data.intersected_gm_before;
            intersected_gm_after{i}(:,j) = repos_data.intersected_gm_after;
            shift_distances{i}(:,j) = repos_data.shift_distances;  
            intersect_lower_boundary_before{i}(:,j) = repos_data.intersect_lower_boundary_before;
            intersect_lower_boundary_after{i}(:,j) = repos_data.intersect_lower_boundary_after;
            inside_gm_before{i}(:,j) = repos_data.inside_gm_before;
            inside_gm_after{i}(:,j) = repos_data.inside_gm_after;
            
            d2layer_before{i}(:,j) = repos_data.d2layer_before;
            repositioned_cells(cell_ids_lin == cell_ij) = 1; % data for this cell found
            
        else
           fprintf('Repositioned data file for %s (cell %g) does not exist\n',...
                    cell_model_names_ij,cell_ij); 
        end
   end
end
fprintf('Reposition data found for %g of %g cells\n',sum(repositioned_cells),length(cell_ids_lin)); 
fprintf('Cell data missing for cell_ids:\n'); 
fprintf('%g\n',find(cell_ids_lin(~repositioned_cells))); 
% Add to NeuronPop
NeuronPop.cell_origins = cell_origins; 
NeuronPop.cell_normals = cell_normals; 
NeuronPop.phis = phis; 
NeuronPop.intersected_gm_before = intersected_gm_before; % remaining intersecting neurons 
NeuronPop.intersected_gm_after = intersected_gm_after; % original intersecting neurons
NeuronPop.shift_distances = shift_distances; 
NeuronPop.repositioned = 1; 
NeuronPop.mode = mode;
NeuronPop.inside_gm_before = inside_gm_before;
NeuronPop.inside_gm_after = inside_gm_after;
NeuronPop.stay_in_layer = stay_in_layer;
NeuronPop.intersect_lower_boundary_before = intersect_lower_boundary_before;
NeuronPop.intersect_lower_boundary_after = intersect_lower_boundary_after;
NeuronPop.intersect_upper_boundary_after = intersect_upper_boundary_after;
NeuronPop.intersect_upper_boundary_before = intersect_upper_boundary_before;
NeuronPop.d2layer_before = d2layer_before;
if in.over_write   
   nrn_pop_file = fullfile(nrn_pop_dir,[nrn_pop_name '_' nrn_model_ver '.mat']); 
   fprintf('Over-writing existing NeuronPop file %s\n',nrn_pop_file); 
else    
    nrn_pop_file = fullfile(nrn_pop_dir,[nrn_pop_name '_' nrn_model_ver '_' mode '_r.mat']); % add r to indicate repositioned population data
    fprintf('Saving NeuronPop %s\n',nrn_pop_file); 
end
save(nrn_pop_file,'NeuronPop'); 
end