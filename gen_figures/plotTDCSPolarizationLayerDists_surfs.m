%% Plot polarization distributions on layer surfaces
% Default settings generate all panels in Fig. 3
Efield_names = {'M1-SO_7x5_rect','M1-HD_4x1_3cm'}; % Name of Efield solutions to plot
amps = [1.8,2.0]; % current amplitude for each Efield
comp_types = {'axon','basal','soma','apical'}; % compartment types
clims_all = {[-0.36 0.45];[-0.36 0.45];[-0.36 0.45];[-0.36 0.45]}; % color axes limits
% clims_all = {[];[];[];[]}; % leave empty to set automatically
save_figs = 0; % set to 1 to save
% By default colorbars not included, add using:
% >> colorbar 

%%
layer_set_num = 10; 
[mat_dir,data_dir] = addPaths_tDCSsim;
params = struct(); 
params.model_prefix_params.layer_set_num = layer_set_num;
params.model_prefix_params.nrn_model_ver = 'maxH';
params.model_prefix_params.dur = 500; 
params.nrn_pop_names = 'nrn_pop1';
args = struct(); 
args.data_fold = fullfile(data_dir,'nrn_sim_data'); 
args.comp_func = 'maxabs';  % function to apply within compartments, see apply_func for options
args.cell_func = 'median';  % function to apply between clones at each location
args.z_lims = [60 90.85]; 
args.shift_dir = [-50,0,0]; 
args.save_fig = save_figs; 
args.fig_fold = fullfile(mat_dir,'figures');
if isnumeric(args.comp_func)
    comp_func_str = num2str(args.comp_func,'%.3f_quant');
else
    comp_func_str = args.comp_func;
end

%% Generate figures
med_pols_all = cell(length(Efield_names),length(comp_types)); % median polarization distributions
for i = 2:length(Efield_names)            
    model_prefix_pre = sprintf('tdcs_r4_v3_%.1f',amps(i));
    params.model_prefix_params.Efield_name = Efield_names{i};     
    params.model_prefix_params.amp = amps(i); 
    params.name_style = 2;
    model_prefix = getNrnSimDataFileName(model_prefix_pre,params);    
    data_struct = load(fullfile(args.data_fold,[model_prefix '.mat']));
    args.deltaVms = data_struct.deltaVms;
    args.cell_model_names = data_struct.cell_model_names;
   for j = 1:length(comp_types)       
       args.comp_type = comp_types{j};  
       if strcmp(args.comp_type,'soma')
           args.sectype = 'all';
       else
           args.sectype = 'terminal';
       end
       args.clims = clims_all{j};
       args.fig_name = sprintf('pol_sameC_ls%g_%s_%s_%s_%s',layer_set_num,...
           args.comp_type,args.sectype,comp_func_str,model_prefix);
       med_pol = plotTDCSpolLayersFig(layer_set_num,nrn_model_ver,model_prefix,args);                      
       med_pols_all{i,j} = med_pol;
        if isempty(clims_all{j})
            clims_all{j} = clim; 
        end
   end
end