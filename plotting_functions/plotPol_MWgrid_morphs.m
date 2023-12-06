function plotPol_MWgrid_morphs(save_fig,nrn_model_ver,comp_func,comp_type,...
                                sectype,plot_comp_type,varargin)
% plotPol_MWgrid_morphs Plots grid of threshold-direction maps (Mollweide projection) 
% Data should already be generated and saved in
% nrn_sim_data/<model_prefix>.mat
if nargin == 0
    save_fig = 0;
    nrn_model_ver = 'maxH';
    comp_func = 'max'; % apply within cell (across compartments), e.g. 'max' or 'maxabs'
    comp_type = 'soma'; % e.g. 'soma', 'axon', or 'dendrite'
    sectype = 'all'; % e.g. 'all' or 'terminal'
    plot_comp_type = 'all'; % same options as comp_type (soma_dendrite or soma_axon only working options)
end
[mat_dir,data_dir] = addPaths_tDCSsim; 

dur = 500; % ms
dtheta = 15; % deg
dphi = 15; % deg
cell_ids = 1:25;
% model_prefix = sprintf('utms_%s_w%g_dth5_dph5',nrn_model_ver,tms_mode); 
cell_lw = 1.5; 
% Plot settings
in.model_prefix_pre = 'upol'; 
in.order_by_maxpol = 1; 
in.mw_settings.display_scale = 'lin';
in.mw_settings.colorbar_auto=1;
in.mw_settings.colorbar_norm=0;
in.mw_settings.colorbar_min=-0.2;
in.mw_settings.colorbar_max=0.2;
in.mw_settings.colorbar_on=0; 
in.mw_settings.title_on=0; 
in.mw_settings.display_grid=0; 
in.mw_settings.plot_point = 'max';
in = sl.in.processVarargin(in,varargin); 
mw_settings = in.mw_settings;
order_by_maxpol = in.order_by_maxpol;
%% Load data
data_fold = fullfile(data_dir,'nrn_sim_data'); % save .mat file here
% Load composite data file (rax=0 for all cells, except rax=2 for cells 32-36,52-56)
model_prefix = sprintf('%s_%s_%gms_dth%g_dph%g',in.model_prefix_pre,nrn_model_ver,dur,dtheta,dphi); 
data = load(fullfile(data_fold,[model_prefix '.mat']));
thetas = data.thetas{1};
phis = data.phis{1}; 
deltaVms = data.deltaVms;
deltaVms1 = processNrnCompData(deltaVms,cell_ids,nrn_model_ver,comp_func,comp_type,sectype); 
%% reorder
if length(order_by_maxpol) == 1 && order_by_maxpol == 1
    deltaVms1(cellfun(@isempty,deltaVms1)) = {nan}; % replace with nan for reshaping below
    max_deltaVms = reshape(cellfun(@(x) max(abs(x)),deltaVms1),5,5); 
    reorder = reshape(1:25,5,5); 
    for i = 1:5
       [~,indi] = sort(max_deltaVms(:,i),'descend');
       reorder(:,i) = reorder(indi,i); 
    end
    % reorder = [1:5:25,2:5:25,3:5:25,4:5:25,5:5:25];
    reorder = reorder'; reorder = reorder(:); 
    cell_ids = cell_ids(reorder); 
    % cell_model_names = cell_model_names(reorder); 
    deltaVms1 = deltaVms1(reorder); 
    deltaVms = deltaVms(reorder); 
elseif length(order_by_maxpol) == length(cell_ids) % input vector for reordering    
    cell_ids = cell_ids(order_by_maxpol); 
    % cell_model_names = cell_model_names(reorder); 
    deltaVms1 = deltaVms1(order_by_maxpol); 
    deltaVms = deltaVms(order_by_maxpol); 
end
if strncmp('max',comp_func,3)
    plot_ind = 'max';
else
    plot_ind = 'min';
end
%% Plot
fig = plotMWmorphGrid(deltaVms1,thetas,phis,cell(size(deltaVms1)),...
    cell_ids,nrn_model_ver,mw_settings,'plot_ind',plot_ind,'comp_data',deltaVms,...
    'lw',cell_lw,'plot_comp_type',plot_comp_type);
for i = (length(fig.Children)/2 + 1):length(fig.Children) % just MW maps
    colormap(fig.Children(i),bluewhitered(1000,fig.Children(i)));    
   if ~mw_settings.colorbar_auto
       if strcmp(mw_settings.display_scale,'log')
            caxis(fig.Children(i),log10([mw_settings.colorbar_min mw_settings.colorbar_max])) % log scale
       else
           caxis(fig.Children(i),[mw_settings.colorbar_min mw_settings.colorbar_max]) % linear scaling
       end
   end
end
% adjust cell 25 plot
if strcmp(nrn_model_ver,'maxH')
    fig.Children(2).XLim = [-1000 1400];
end
% adjust cell 8 plot
% fig.Color = 0.4*ones(1,3); 
fig.Color = 'w';
fig.Position = [1 1/6 0.75 0.7444]; 
% fig.Units = 'inches';
% fig.Position(3:4) = [20.0000   90.0972];
% save fig
if save_fig
    if length(order_by_maxpol) == 1
       fig_name = sprintf('%s_mwgrid_%g_cauto%g_%s_%s_%s',model_prefix,...
           order_by_maxpol,mw_settings.colorbar_auto,comp_type,sectype,comp_func);
    else
        fig_name = sprintf('%s_mwgrid_cord_cauto%g_%s_%s_%s',model_prefix,...
                        mw_settings.colorbar_auto,comp_type,sectype,comp_func);
    end
    printFig(fig,fullfile(mat_dir,'figures'),fig_name,'formats',{'png'},...
                'resolutions',{'-r350'});
end