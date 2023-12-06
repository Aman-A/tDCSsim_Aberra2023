function Maps = getMapsCompsCellM(model_prefix,layer_set_num,Efield_name,nrn_pop_names,...
                                nrn_model_ver,map_type,cell_ids,comp_type,sectype,func,varargin)
%GETMAPS Use uniform E threshold-direction or polarization-direction map
%and Efield distribution to to estimate thresholds or polarizations of
%neurons, respectively
%Uses E at cell body
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
in.load_maps = 1; % default load if exists
in.save_maps = 1; % default save after generating
in.reposition_mode = 'off';
in = sl.in.processVarargin(in,varargin);
% Set paths
[mat_dir,data_dir] = addPaths_tDCSsim;
map_fold = fullfile(data_dir,'nrn_sim_data','map_data');
% Get Map file name
if iscell(nrn_pop_names)
    num_pops = length(nrn_pop_names);
elseif ischar(nrn_pop_names)
    nrn_pop_names = {nrn_pop_names};
    num_pops = 1;
end
if isnumeric(func)
    func_str = num2str(func,'%.3f_quant');
else
    func_str = func;
end
if num_pops == 1
    map_filename = sprintf('%s_ls_%g_E_%s_P_%s_cell_%s_%s_%s',model_prefix,layer_set_num,...
        Efield_name,nrn_pop_names{1},...
        func_str,comp_type,sectype);
else
    map_filename = sprintf('%s_ls_%g_E_%s_P_%s-%s_cell_%s_%s_%s',model_prefix,layer_set_num,...
        Efield_name,nrn_pop_names{1},nrn_pop_names{end},...
        func_str,comp_type,sectype);
end
map_file_path = fullfile(map_fold,[map_filename '.mat']);
% Load file if exists
if in.load_maps && exist(map_file_path,'file')
    fprintf('Loading maps from %s\n',map_file_path); 
    Maps_data = load(map_file_path); 
    Maps = Maps_data.Maps;
    return; 
end
% Otherwise proceed to make Maps
nrn_phis_all = cell(num_pops,1);
layersE_all = cell(num_pops,1);
for pop = 1:num_pops
    NeuronPop = loadNeuronPop(layer_set_num,nrn_pop_names{pop},nrn_model_ver);
    nrn_phis_all{pop} = NeuronPop.phis;
    if ~strcmp(in.reposition_mode,'off') && ~strcmp(in.reposition_mode,'none')
        layersE = loadLayers(layer_set_num,'opt','layersE','Efield_name',Efield_name,...
                        'reposition_mode',in.reposition_mode,'nrn_model_ver',nrn_model_ver,...
                        'nrn_pop_name',nrn_pop_names{pop});
        layersE_all{pop} = layersE; % load unique layersE for each population
    else
        if pop == 1 % load layersE just for first population, all are identical
            layersE = loadLayers(layer_set_num,'opt','layersE','Efield_name',Efield_name); 
            layersE_all{pop} = layersE; 
        else
            layersE_all{pop} = layersE_all{1}; 
        end
    end
end
fprintf('Loaded phis from %g NeuronPops\n',num_pops);
switch map_type
    case 'polarization'
        % get indices to include in loaded data
        % extract relevant data
        cell_ids_lin = [cell_ids{:}]; % turn cell array into linear vector
        % Load data
        deltaVmMap_data = load(fullfile(mat_dir,'nrn_sim_data',[model_prefix '.mat']));
        cell_model_names = deltaVmMap_data.cell_model_names;
        cell_ids_lin_sim = getCellIds(cell_model_names); % cell_ids in linear vector
        [~,inds_include,~] = intersect(cell_ids_lin_sim,cell_ids_lin);
        cell_model_names = cell_model_names(inds_include);
        deltaVmsU = deltaVmMap_data.deltaVms(inds_include);
        thetasU = deltaVmMap_data.thetas(inds_include);
        phisU = deltaVmMap_data.phis(inds_include);

        % extract polarization data for each direction as func applied to
        % comp_type/sectype compartments

        deltaVmsUcompStat  = processNrnCompData(deltaVmsU,cell_ids_lin,...
                                                nrn_model_ver,func,comp_type,sectype);
        % Make interpolant maps consisting of func applied to comp_type polarization
        % values at each direction (theta/phi) uniform E -field
        deltaVms_int_grids = makeMWgriddedInterp(thetasU,phisU,deltaVmsUcompStat);
        num_cells = length(cell_model_names);
        deltaVms = cell(1,num_cells); % polarization data at each position in layers
        num_layers =  length(cell_ids);
        E_vals = cell(num_layers,1);        
        fprintf('Interpolating within map %s polarizations in comp_type: %s, sectype: %s, using E at somas\n',...
            func_str,comp_type,sectype);
        % Use map to get deltaVms in layer for each cell/rotation
        for li = 1:num_layers
            num_cells_in_layer = length(cell_ids{li});                        
            E_vals{li} = cell(1,num_cells_in_layer);
            for cj = 1:num_cells_in_layer
                cell_model_name = cellModelNames(cell_ids{li}(cj));
                [~,~,cell_ind] = intersect(cell_model_name,cell_model_names);                                
                deltaVms{cell_ind} = zeros(size(layersE(li).Efield,1),num_pops);
                % transform E to local neuron space (s-d axis aligned to [0 0 1])
                E_vals{li}{cj} = cell(1,num_pops);
                for pop = 1:num_pops
                    if size(layersE_all{pop}(li).Efield,3) == num_cells_in_layer
                        Eicp = layersE_all{pop}(li).Efield(:,:,cj);   
                        cell_normals = layersE_all{pop}(li).cell_normals(:,:,cj);
                    else
                        Eicp = layersE_all{pop}(li).Efield;
                        cell_normals = layersE_all{pop}(li).cell_normals;
                    end
                    Ei_loc = zeros(size(Eicp));
                    for posk = 1:size(Eicp,1)
                        Ei_loc(posk,:) = reorientEfield(cell_normals(posk,:),nrn_phis_all{pop}{li}{cj}(posk),Eicp(posk,:));
                    end
                    E_vals{li}{cj}{pop} = Ei_loc;
                    [phis,lambdas,Emag_fem] = cart2sph(Ei_loc(:,1),Ei_loc(:,2),Ei_loc(:,3)); % convert to spherical coords
                    Etheta_fem = 90-lambdas*180/pi; % convert to polar angle (angle from z) in deg
                    Ephi_fem = phis*180/pi; % convert to deg
                    Ephi_fem(Ephi_fem >= 180) = Ephi_fem(Ephi_fem >= 180) - 360; % set to -180 to 180 deg
                    deltaVm_int_grid_ij = deltaVms_int_grids{cell_ind};
                    deltaVms{cell_ind}(:,pop) = Emag_fem.*deltaVm_int_grid_ij(Etheta_fem,Ephi_fem);
                end
            end
        end
        fprintf('Finished calculating deltaVms in layers\n')
        if length(layersE_all) == 1
            layersE_all = layersE_all{1}; 
        end
        Maps.model_prefix = model_prefix;
        Maps.nrn_model_ver = nrn_model_ver;
        Maps.layersE = layersE_all;
        Maps.E_loc = E_vals; % E-fields for each cell model and rotation
        Maps.Efield_name = layersE(1).Efield_name;
        Maps.nrn_pop_names = nrn_pop_names;
        Maps.nrn_phis_all = nrn_phis_all;
        Maps.cell_ids = cell_ids;
        Maps.deltaVmsUcompStat = deltaVmsUcompStat; % polarization stat (func) in comp_type/sectype
        Maps.deltaVm_int_grids = deltaVms_int_grids; % interpolation grids for deltaVmsUcompStat
        Maps.deltaVms = deltaVms; % polarization stat in FEM E-field at each location
        Maps.cell_model_names = cell_model_names;
        Maps.func = func; % function applied to deltaVmsU
        Maps.comp_type = comp_type; % compartment type
        Maps.sectype = sectype;  % section type
        %% Get deltaVms after taking weighted mean of statistic across clones and rotations
        [~,~,layer_mtype_names,etypes] = cellfun(@(x) cellNameParser(x),cell_model_names,'UniformOutput',0);
        layer_metype_names = join([layer_mtype_names,etypes],'_');
        [layer_metype_namesU,~,indsb] = unique(layer_metype_names);
        num_unique_metypes = length(layer_metype_namesU); % number of maps after averaging across clones/rotations
        %         w = cell(size(cell_model_names)); % same structure as cell_ids, allows for variable number of clones per layer
        num_comps = zeros(num_cells,1);
        for i = 1:num_cells
            celli = loadCellData(cell_ids_lin(i),nrn_model_ver,{'comp_types','sectypes','secnames'});
            ind = getCellCompInds(celli,comp_type,sectype);
            num_comps(i) = sum(ind);
        end
        % Compute weighted mean of deltaVm to uniform field
        % makes Map struct data layer specific
        deltaVms_weighted = cell(1,num_unique_metypes); % one vector per layer, combining data across clones/rotations
        cell_model_names_weighted = cell(num_unique_metypes,1); % full name of 1st clone in each set, using full name convention (for use with other functions)
        for i = 1:num_unique_metypes
            clone_inds = find(indsb == i); % indices in full vector for these clones of same me-type
            cell_model_names_weighted{i} = cell_model_names{clone_inds(1)};
            deltaVms_clonesi = deltaVms(clone_inds); % polarization data for clones (num_positions x num_pops)
            if strcmp(comp_type,'soma') % exception for soma case (single comp per cell/rotation)
                % apply function across clones/rotations by separating each
                deltaVms_weighted{i} = apply_func(cell2mat(deltaVms_clonesi), func);
                if i == 1
                    fprintf('Applying func: %s to somatic polarizations across %g clones and %g rotations\n',...
                            func_str, length(clone_inds), num_pops);
                end
            else
                % weights are number of compartments divided by sum of all
                % compartments and number of rotations
                weights = num_comps(clone_inds)/(sum(num_comps(clone_inds))*num_pops);
                weights_cell = mat2cell(weights,ones(length(weights),1))';
                % sum of weighted polarization values across all clones/rotations
                % scale each column (pop/rotation) by same weight within cell
                deltaVms_clonesi_weighted = cellfun(@(x,y) x*y,...
                    deltaVms_clonesi,weights_cell,'UniformOutput',0);
                deltaVms_weighted{i} = sum(cell2mat(deltaVms_clonesi_weighted),2); %
            end
        end
        % Add weighted mean to Maps struct - remove empty elements
        Maps.deltaVms_weighted = deltaVms_weighted; % weighted means in each element across clones/rotations after applying func within compartments
        Maps.cell_model_names_layer_weighted = layer_metype_namesU; % name of m/e type (cell type) (no clone number)
        Maps.cell_model_names_weighted = cell_model_names_weighted; % name of first clone for each set of weighted means across clones & rotations    
end
if in.save_maps    
    if ~exist(map_fold,'dir')
        mkdir(map_fold)
    end
    save(map_file_path,'Maps');
    fprintf('Saved Maps to %s\n',map_filename);
end