function deltaVm_roi = plotPolInROIViolin(layer_set_num,nrn_model_ver,model_prefixes,varargin)
%PLOTPOLINROIVIOLIN plot violin plots of threshold distributions within FDI
%representation with experimental data overlaid for Fig. 6
%
% AUTHOR    : Aman Aberra
mat_dir = addPaths_tDCSsim;
if nargin == 0
    layer_set_num = 1;
    nrn_model_ver = 'maxH';
    % model name
    model_prefix_pre = 'tdcs';
    dur = 500;
    Efield_name = 'M1-SO_conv_1mA';
    nrn_pop = 'nrn_pop1';
    model_prefixes = {sprintf('%s_%s_%gms_ls_%g_E_%s_P_%s',model_prefix_pre,nrn_model_ver,dur,...
                                layer_set_num,Efield_name,nrn_pop)};
end
in.data_fold = fullfile(mat_dir,'nrn_sim_data');
in.cell_ids = {1:5;6:10;11:15;16:20;21:25};
in.model_names = {'HD'};
in.plot_region_name = ''; % 'FDI_rep_inds', 'handknob_inds'
in.comp_func = 'max'; % apply within cell (across compartments), e.g. 'max' or 'maxabs'
in.comp_type = 'soma'; % e.g. 'axon' or 'dendrite'
in.sectype = 'all'; % e.g. 'all' or 'terminal'
in.deltaVms_all_models = {};
in.cell_model_names = {};
in.num_pops = 1;
% Plot settings
in.font_size = 12;
in.font_name = 'Times';
in.yaxlim = []; % A/us
in.y_tick_marks = [];
in.yscale = 'linear';
in.colors = jet(length(in.cell_ids));
in.colors(4,:) = [0.5961 0.3059 0.6392];
in.fig_units = 'centimeters';
in.fig_size = [14.9,11.2];
in.y_label = 'Polarization (mV)';
in.violin_lw = 2; % line width of mean/median line in violin plot
in.plot_pts = 0; 
in.save_fig = 0;
in.fig_fold = fullfile(mat_dir,'figures');
in.fig_name = ['polviol_' model_prefixes{1}];
in = sl.in.processVarargin(in,varargin);
%% Load data
% Polarization data
% Load threshold data
num_models = length(model_prefixes);
if isempty(in.deltaVms_all_models)
    fprintf('Loading data...\n')
    deltaVms_all_models = cell(num_models,1);
    % loads data
    for i = 1:num_models
        data_structi = load(fullfile(in.data_fold,[model_prefixes{i} '.mat']));
        deltaVms_all_models{i} = data_structi.deltaVms;
    end
    cell_model_names = data_structi.cell_model_names;
else
    fprintf('Data already loaded\n')
    deltaVms_all_models = in.deltaVms_all_models;
    cell_model_names = in.cell_model_names;
end
%% Extract data in desired cell models
plot_cell_ids = [in.cell_ids{:}];
[~,~,data_inds] = intersect(cellModelNames(plot_cell_ids),cell_model_names);
for i = 1:num_models
    deltaVms_all_models{i} = deltaVms_all_models{i}(data_inds);
end
cell_model_names = cell_model_names(data_inds);
plot_layers = cellfun(@(x) ~isempty(x),in.cell_ids,'UniformOutput',1);
%% Apply function in desired compartments
if in.num_pops == 1 % for num_pop > 1, data should already have single value per population
    % preload cell data 
    cell_data_all = cell(length(plot_cell_ids),1); 
    for ii = 1:length(plot_cell_ids)
        cell_data_all{ii} = loadCellData(plot_cell_ids(ii),nrn_model_ver,...
                                        {'comp_types','sectypes','secnames'});
                                        
    end    
    fprintf('Processing data...\n');
    for i = 1:num_models
        % Process polarization values to get single value for each neuron
        deltaVms_all_models{i} = processNrnCompData(deltaVms_all_models{i},...
            plot_cell_ids,nrn_model_ver,in.comp_func,in.comp_type,in.sectype,...
            'cell_data',cell_data_all);
    end
    fprintf('Done\n'); 
end
%% Load layerROI
layers = loadLayers(layer_set_num);
if isempty(in.plot_region_name) || strcmp(in.plot_region_name,'all')
    % use data from all neurons    
    num_elems = [layers.num_elem];
    ROIi = arrayfun(@(x) true(x,1),num_elems,'UniformOutput',0)';
else % load ROI for extracting data in plot_region_name
    if isfield(layers(1),'region_labels')
       if any(strcmp(in.plot_region_name,layers(1).region_labels))
           % check if input plot_region_name is one of the region_labels,
           % otherwise assume it's a manually defined layerROI
            use_label_ind = find(strcmp(in.plot_region_name,layers(1).region_labels));       
       else
           error('%s is not defined a region label',in.plot_region_name)
       end
    end
    ROIi = cell(length(layers),1);
    for i = 1:length(layers)
        ROIi{i} = layers(i).region_finds == use_label_ind;
    end
   
end
%% Plot violin plots
fig = figure('Color','w');
deltaVm_roi = plotDistsModelLayerViolin(deltaVms_all_models,ROIi,...
                            in.cell_ids,cell_model_names,'colors',in.colors,...
                            'violin_lw',in.violin_lw,'plot_pts',in.plot_pts);
ax = gca;
ax.YScale = in.yscale;
if ~isempty(in.yaxlim)
    ax.YLim = in.yaxlim;
end
if ~isempty(in.y_tick_marks)
    ax.YTick = in.y_tick_marks;
end
box(ax,'off');
if num_models > 1
    delete(fig.Children(1)); % remove legend
end
ax.FontName = in.font_name;
ax.FontSize = in.font_size;
if num_models > 1
    ax.XTickLabel = in.model_names;
else
   layer_names = {'L1','L2/3','L4','L5','L6'};
   ax.XTickLabel = layer_names(plot_layers);
end
fig.Units = in.fig_units;
fig.Position(3:4) =  in.fig_size;
ax.YGrid = 'on'; ax.YMinorGrid = 'off';
ax.XColor = 'k'; ax.YColor = 'k';
ylabel(in.y_label);
if in.save_fig
    printFig(fig,in.fig_fold,in.fig_name,'formats',{'fig','png'},...
                'resolutions',{[],'-r600'});   
end
end
