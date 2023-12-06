function layers_or_layersE = loadLayers(layer_set_num,varargin)
% LOADLAYERS Load layers struct or layersE struct
%
if nargin == 0
   layer_set_num = 1;
end
mat_dir = addPaths_tDCSsim;
in.opt = 'layers';
in.Efield_name = '';
in.reposition_mode = 'off';
in.nrn_pop_name = ''; % required if repositioning was used
in.nrn_model_ver = ''; % required if repositioning was used
in = sl.in.processVarargin(in,varargin);
if nargin < 1
   layer_set_num = 1;
end
layer_set_name = sprintf('layer_set_%g',layer_set_num);
if strcmp(in.opt,'layers')
    layer_file_name = sprintf('%s.mat',layer_set_name);
    is_layers = true;
elseif strcmp(in.opt,'layersP') || strcmp(in.opt,'layersp')
    % layersP includes with both cell placement layer and boundaries
    layer_file_name = sprintf('%sp.mat',layer_set_name);
    is_layers = true;
elseif strcmp(in.opt,'layersE') || strcmp(in.opt,'layerse')
    % layersE includes E-field data for Efield_name within layers
    if ~strcmp(in.reposition_mode,'off') && ~strcmp(in.reposition_mode,'none')
        layer_file_name = sprintf('%s_E_%s_%s_%s_%s.mat',layer_set_name,in.Efield_name,...
            in.reposition_mode,in.nrn_model_ver,in.nrn_pop_name);
    else
        layer_file_name = sprintf('%s_E_%s.mat',layer_set_name,in.Efield_name);
    end
    is_layers = false; % layersE
else
    error('%s not a valid option',in.opt);
end

%% Load
% fprintf('Loading %s\n',layer_file_name);
layer_table = readtable(fullfile(mat_dir,'layer_sets.csv'));
mesh_name = layer_table.mesh_name{layer_table.ID==layer_set_num};
roi_name = layer_table.roi_name{layer_table.ID==layer_set_num};
mesh_roi_name = [mesh_name '_' roi_name];
mesh_roi_folder = fullfile(mat_dir,'output_data','layer_data',mesh_name,mesh_roi_name);
layer_data = load(fullfile(mesh_roi_folder,layer_set_name,layer_file_name));

if is_layers
    layers_or_layersE = layer_data.layers;
else
    layers_or_layersE = layer_data.layersE;
end