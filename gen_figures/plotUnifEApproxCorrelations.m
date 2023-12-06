%% ** Plot Fig 8 and Fig S10 **
% Plots scatter plots showing correlation between simulated polarizations and
% estimated polarizations with uniform E-field approximation (Fig S10)
% and bar plots of error statistics (Fig 8)
layer_set_num = 10;
reposition_mode = 'both';
model_prefix_pres = {'tdcs_r4_v3','tdcs_r4_v3'};
nrn_model_ver = 'maxH';
dur = 500;     
amps = [1.8,2.0]; 
Efield_names = {'M1-SO_7x5_rect','M1-HD_4x1_3cm'};
nrn_pop_names = {'nrn_pop1'}; 
comp_types = {'axon','basal','soma','apical'};
sectype = 'all';
sectypes = {'terminal','terminal','all','terminal'};
comp_func = 'maxabs'; % do processing externally
params = struct(); 
params.model_prefix_params.layer_set_num = layer_set_num;
params.model_prefix_params.nrn_model_ver = 'maxH';
params.model_prefix_params.dur = 500; 
params.nrn_pop_names = nrn_pop_names;
[mat_dir,data_dir] = addPaths_tDCSsim;

cell_ids = {1:5;6:10;11:15;16:20;21:25};
data_fold = fullfile(data_dir,'nrn_sim_data'); 
cell_func = 'median';  
% Map settings
map_model_prefix = 'upol_maxH_500ms_dth15_dph15';
map_Efield_name = 'M1-SO_7x5_rect';
map_nrn_pop_names = nrn_pop_names; 
map_type = 'polarization'; 
apic_cell_ids = {[],6:10,[],16:20,21:25}; % no interneurons
apic_cell_ids_lin = [apic_cell_ids{:}];
cell_ids_lin = [cell_ids{:}]; 
cell_model_names = cellModelNames(cell_ids_lin);
%% Load neuron simulation data and Map data, process for compartment stats
deltaVms_all = cell(length(Efield_names),length(comp_types)); % nrn simulation data
deltaVmsM_all = cell(length(Efield_names),length(comp_types)); % map estimation data
if isnumeric(comp_func)
    comp_func_str = num2str(comp_func,'%.3f_quant');
else
    comp_func_str = comp_func;
end
processed_data_fold = fullfile(mat_dir,'gen_figures','processed_data');
if ~exist(processed_data_fold,'dir')
    mkdir(processed_data_fold)
end
processed_sim_data_file = fullfile(processed_data_fold,...
                                  sprintf('ls%g_%s_%s_%gE.mat',layer_set_num,...
                                  comp_func,sectype,length(Efield_names))); 
if ~exist(processed_sim_data_file,'file')
    for i = 1:length(Efield_names)
        Efield_name = Efield_names{i}; 
        model_prefix_pre = sprintf('tdcs_r4_v3_%.1f',amps(i));
        params.model_prefix_params.Efield_name = Efield_name;     
        params.model_prefix_params.amp = amps(i); 
        params.name_style = 2; % 1 - old simulations
        model_prefix = getNrnSimDataFileName(model_prefix_pre,params);    
        data_struct = load(fullfile(data_fold,[model_prefix '.mat']));
        fprintf('Loaded %s\n',model_prefix); 
        % layersE for map estimation with this E-field
        layersE = loadLayers(layer_set_num,'opt','layersE','Efield_name',Efield_name,...
                            'reposition_mode',reposition_mode,'nrn_model_ver',nrn_model_ver,...
                            'nrn_pop_name',nrn_pop_names{1}); 
        for j = 1:length(comp_types)
            deltaVmsij = processNrnCompData(data_struct.deltaVms,cell_ids_lin,...
                                            nrn_model_ver,comp_func,comp_types{j},...
                                            sectypes{j},'print_level',1); 
            deltaVms_all{i,j} = deltaVmsij;    
            % Load map and generate estimation data
            if strcmp(comp_types{j},'apical')
                cell_idsj = apic_cell_ids;
            else
                cell_idsj = cell_ids;
            end
            % Get Maps for this compartment type and E-field
            Maps = getMapsCompsCellM(map_model_prefix,layer_set_num,map_Efield_name,...
                                    map_nrn_pop_names,nrn_model_ver,map_type,...
                                    cell_idsj,comp_types{j},sectype,comp_func,...
                                    'reposition_mode',reposition_mode,'load_maps',0);
            deltaVmsMij = mapEstLayers(Maps,layersE,cell_idsj,...
                                    nrn_pop_names,'scale_E',amps(i));
            deltaVmsM_all{i,j} = deltaVmsMij;
            
        end
    end
    save(processed_sim_data_file,'deltaVms_all','deltaVmsM_all'); 
    fprintf('Saved processed simulation and Map data to %s\n',processed_sim_data_file); 
else
    proc_data = load(processed_sim_data_file); 
    deltaVms_all = proc_data.deltaVms_all; 
    deltaVmsM_all = proc_data.deltaVmsM_all;
    fprintf('Loaded processed simulation and Map data from %s\n',processed_sim_data_file); 
end
%% Collapse across layers
fprintf('Collapsing across layers\n'); 
deltaVms_layer_all = cell(length(Efield_names),length(comp_types)); % 
deltaVmsM_layer_all = cell(length(Efield_names),length(comp_types)); % 
for i = 1:length(Efield_names)
   for j = 1:length(comp_types)       
       deltaVmsij = deltaVms_all{i,j}; 
       if strcmp(comp_types{j},'apical')
           deltaVmsMij = cell(1,length(cell_ids_lin)); 
           deltaVmsMij(apic_cell_ids_lin) = deltaVmsM_all{i,j};            
       else
           deltaVmsMij = deltaVmsM_all{i,j}; 
       end
       deltaVms_layer_all{i,j} =  calcDataLayers(deltaVmsij,cell_model_names,...
                                                cell_ids,'func','none');
       deltaVmsM_layer_all{i,j} =  calcDataLayers(deltaVmsMij,cell_model_names,...
                                                   cell_ids,'func','none');       
   end
end
%% Plot scatter plots
% colors = jet(length(cell_ids));
% colors(4,:) = [0.5961 0.3059 0.6392];
sc_opts = struct;
sc_opts.median_in_layer = 0;
sc_opts.save_fig = 0; 
sc_opts.fig_size = [4.403 7.5556];
sc_opts.mcol = lines(5);
sc_opts.fig_fold = fullfile(mat_dir,'figures'); 
rsq_all = cell(length(Efield_names),length(comp_types)); 
p_all = cell(length(Efield_names),length(comp_types)); 
med_ape_all = cell(length(Efield_names),length(comp_types)); 
mane2_all = cell(length(Efield_names),length(comp_types));  % mean abs norm error by 0.975 quantile
mape_all = cell(length(Efield_names),length(comp_types));  % mean abs norm error by 0.975 quantile
for i = 1:length(Efield_names)
    Efield_name = Efield_names{i};
    fprintf('E: %s\n',Efield_name); 
    for j = 1:length(comp_types)
        comp_type = comp_types{j};
        fprintf(' comp_type: %s\n',comp_type); 
        x_all = deltaVms_layer_all{i,j};
        y_all = deltaVmsM_layer_all{i,j};
        sc_opts.fig_name = sprintf('unifECorrel_ls%g_%s_%s_%s',...
                                    layer_set_num,comp_type,comp_func_str,...
                                    Efield_name);
        [rsq_all{i,j},p_all{i,j},med_ape_all{i,j},mane2_all{i,j},mape_all{i,j}] = ...
            plotScatterPlotGrid(x_all,y_all,sc_opts); 
    end
end
%% Plot bar plots 
save_rsq_fig = 0; 
alpha = 0.005; 
bar_opts = struct; 
bar_opts.cell_labels = {'L1','L2/3','L4','L5','L6'}; 
bar_opts.comp_types = comp_types;                   
bar_opts.x_shifts = 0; 
% bar_opts.bar_cols = 0.8*ones(1,3); 
bar_opts.bar_cols = jet(length(cell_ids)); 
bar_opts.bar_cols(4,:) = [0.5961 0.3059 0.6392];
bar_opts.YTick = []; 
bar_opts.bar_width = 0.5; 
bar_opts.y_lim_same = 1;
for i = 1:length(Efield_names)
    Efield_name = Efield_names{i}; 
    med_ape_alli = med_ape_all(i,:); 
    bar_opts.fig_fold = fullfile(mat_dir,'figures');
    bar_opts.fig_name = sprintf('unifEMedAPE_ls%g_%s_%s',...
                                    layer_set_num,comp_func_str,...
                                    Efield_name);
    sig_all = cellfun(@(x) x < alpha,p_all(i,:),'UniformOutput',0);
    sig_all = cell2mat(sig_all); 
    if all(sig_all)
       sig_all = ones(length(med_ape_all{1}),length(comp_types),1);  
    end
    plotRsq_barplot(med_ape_alli,cell(size(med_ape_alli)),cell(size(med_ape_alli)),...
         sig_all,save_rsq_fig,bar_opts);
end
%% Plot MedAPE together
save_rsq_fig = 0; 
% bar_opts.x_shifts = [-2/9,2/9]; 
bar_opts.x_shifts = [-0.18,0.18]; 
bar_opts.bar_cols = [1 0 0;0 0 1]; 
bar_opts.bar_width = bar_opts.x_shifts(2)*2; 
bar_opts.y_lim_same = 1;
bar_opts.YTick = []; 
bar_opts.fig_name = sprintf('unifEMedAPE_ls%g_%s_%s_%s',...
                                    layer_set_num,comp_func_str,...
                                    Efield_names{1},Efield_names{2});
sig_all = ones(length(med_ape_all{1}),length(comp_types),2);  
plotRsq_barplot(med_ape_all(1,:),med_ape_all(2,:),cell(size(med_ape_all(1,:))),...
     sig_all,save_rsq_fig,bar_opts);
%% Plot MAPE together
save_rsq_fig = 0; 
% bar_opts.x_shifts = [-2/9,2/9]; 
bar_opts.x_shifts = [-0.18,0.18]; 
bar_opts.bar_cols = [1 0 0;0 0 1]; 
bar_opts.bar_width = bar_opts.x_shifts(2)*2; 
bar_opts.y_lim_same = 1;
bar_opts.YTick = []; 
bar_opts.fig_name = sprintf('unifEMAPE_ls%g_%s_%s_%s',...
                                    layer_set_num,comp_func_str,...
                                    Efield_names{1},Efield_names{2});
sig_all = ones(length(mape_all{1}),length(comp_types),2);  
plotRsq_barplot(mape_all(1,:),mape_all(2,:),cell(size(mape_all(1,:))),...
     sig_all,save_rsq_fig,bar_opts);
%% Plot MANE together
save_rsq_fig = 0; 
% bar_opts.x_shifts = [-2/9,2/9]; 
bar_opts.x_shifts = [-0.18,0.18]; 
bar_opts.bar_cols = [1 0 0;0 0 1]; 
bar_opts.bar_width = bar_opts.x_shifts(2)*2; 
bar_opts.y_lim_same = 1;
bar_opts.YTick = []; 
bar_opts.fig_name = sprintf('unifEMANE_ls%g_%s_%s_%s',...
                                    layer_set_num,comp_func_str,...
                                    Efield_names{1},Efield_names{2});
sig_all = ones(length(mane2_all{1}),length(comp_types),2);  
plotRsq_barplot(mane2_all(1,:),mane2_all(2,:),cell(size(mane2_all(1,:))),...
     sig_all,save_rsq_fig,bar_opts); 
 