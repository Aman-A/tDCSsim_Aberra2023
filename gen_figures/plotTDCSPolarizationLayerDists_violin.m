%% Plot polarization distributions as violin plots
% Plots Fig. 5 and Fig S7 (region-specific distributions)
comp_types = {'axon','basal','soma','apical'}; % loop over comp_types
lsn = 10; % 3 
nrn_model_ver = 'maxH';
plot_region_names = {'G_precentral','S_central','G_postcentral'}; 
Efield_names = {'M1-SO_7x5_rect','M1-HD_4x1_3cm'};
yaxlims = struct();
yaxlims.soma =  [-0.085 0.085]; % mV
yaxlims.axon = [-1 1]; % mV
yaxlims.basal = [-0.25 0.25]; % mV
yaxlims.apical = [-0.25 0.25]; % mV 
% Constant plot arguments
[mat_dir,data_dir] = addPaths_tDCSsim;
args = struct();
args.save_fig = 0; % set to 1 to save 
args.comp_func = 'maxabs'; % function applied within each cell across compartments
args.sectype = 'all'; % e.g. 'all' or 'terminal'
args.cell_ids = {1:5;6:10;11:15;16:20;21:25};  
args.model_names = {'M1-SO','HD 4x1'}; 
args.font_size = 8;
args.font_name = 'Arial;';
args.y_label = 'Peak polarization (mV)'; 
args.yscale = 'linear';
args.violin_lw = 0.5; 
args.plot_pts = 1; 
args.colors = jet(length(args.cell_ids)); 
args.colors(4,:) = [0.5961 0.3059 0.6392];
args.fig_fold = fullfile(mat_dir,sprintf('figures_%s',nrn_model_ver)); 
args.fig_units = 'inches';
args.fig_size = [2.2 1.8];
args.data_fold = fullfile(data_dir,'nrn_sim_data');
%% Load data with recordings from all comps
amps = [1.8, 2.0];  
model_prefix_pres = arrayfun(@(x) sprintf('tdcs_r4_v3_%.1f',x),amps,'UniformOutput',0);
params = struct(); 
params.model_prefix_params.layer_set_num = lsn;
params.model_prefix_params.nrn_model_ver = nrn_model_ver;
params.model_prefix_params.dur = 500; 
params.nrn_pop_names = 'nrn_pop1'; 
params.name_style = 2;
model_prefixes = cell(length(model_prefix_pres),1); 
for i = 1:length(model_prefix_pres)     
    params.model_prefix_params.Efield_name = Efield_names{i};     
    params.model_prefix_params.amp = amps(i); 
    model_prefixes{i} = getNrnSimDataFileName(model_prefix_pres{i},params);
end
deltaVms_all_models = cell(length(model_prefixes),1);
% loads data
for i = 1:length(model_prefixes)
    data_filei = fullfile(data_dir,'nrn_sim_data',[model_prefixes{i} '.mat']);
    fprintf('Loading %s...\n',data_filei); 
    data_structi = load(data_filei);
    deltaVms_all_models{i} = data_structi.deltaVms;    
end
fprintf('Done loading data\n'); 
cell_model_names = data_structi.cell_model_names;
args.deltaVms_all_models = deltaVms_all_models; 
args.cell_model_names = cell_model_names;
%% Plot peak polarization dist in full ROI (FIG 4)
args.plot_region_name = 'all';
args.fig_size = [2.0 1.8];
pols_all = cell(length(comp_types),1); 
for c = 1:length(comp_types)
    args.comp_type = comp_types{c}; % e.g. 'axon' or 'dendrite'   
    if strcmp(args.comp_type,'apical')
        args.cell_ids = {[];6:10;[];16:20;21:25};  % inhibitory cells don't have apical dendrites 
    else
        args.cell_ids = {1:5;6:10;11:15;16:20;21:25};  
    end
    if c == 1
        args.y_label = 'Peak polarization (mV)';
    else
        args.y_label = '';
    end
    if strcmp(args.comp_type,'soma')        
        args.sectype = 'all';        
    else
        args.sectype = 'terminal';        
    end
    args.yaxlim = yaxlims.(args.comp_type); 
    args.fig_name = sprintf('polViolin_pts%g_ls%g_%s_%s_%s_%s_%s',args.plot_pts,lsn,nrn_model_ver,...
                            args.plot_region_name,args.comp_func,args.comp_type,...
                            args.sectype); 
    pols = plotPolInROIViolin(lsn,nrn_model_ver,model_prefixes,args);
    print_pol_stats(pols,args.comp_func,args.comp_type,args.plot_region_name,...
                    'cutoff',0.5); 
    pols_all{c} = pols;     
end
%% Plot peak polarization dist in each region (Fig S7)
for c = 1:length(comp_types)    
    args.comp_type = comp_types{c}; % e.g. 'axon' or 'dendrite'    
    args.yaxlim = yaxlims.(args.comp_type); 
    if strcmp(args.comp_type,'soma')
        args.sectype = 'all';
    else
        args.sectype = 'terminal';
    end
    if strcmp(args.comp_type,'apical')
        args.cell_ids = {[];6:10;[];16:20;21:25};  % inhibitory cells don't have apical dendrites 
    else
        args.cell_ids = {1:5;6:10;11:15;16:20;21:25};  
    end
    for i = 1:length(plot_region_names)
        if i > 1
           args.y_label = '';  % different regions on same row
        else
            args.y_label = 'Peak polarization (mV)';          
        end
        args.plot_region_name = plot_region_names{i}; 
        fprintf('***plot_region_name = %s:...\n',args.plot_region_name); 
        args.fig_name = sprintf('polViolin_pts%g_ls%g_%s_%s_%s_%s_%s',args.plot_pts,...
                                lsn,nrn_model_ver,args.plot_region_name,...
                                args.comp_func,args.comp_type,...
                                args.sectype); 
        pols = plotPolInROIViolin(lsn,nrn_model_ver,model_prefixes,args);        
        print_pol_stats(pols,args.comp_func,args.comp_type,args.plot_region_name,...
                         'cutoff',0.5); 
    end
end