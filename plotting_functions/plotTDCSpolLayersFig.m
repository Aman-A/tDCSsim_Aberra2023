function pol_data = plotTDCSpolLayersFig(layer_set_num,nrn_model_ver,model_prefix,varargin)
% Plots median polarization across clones and rotations at each location on
% layers 1-5, adjacent to each other, as in Fig. 3c
mat_dir = addPaths_tDCSsim;
if nargin==0
    layer_set_num = 1;
    nrn_model_ver = 'maxH';
    model_prefix_pre = 'tdcs';
    dur = 500;
    Efield_name = 'M1-SO_conv_1mA';
    nrn_pop = 'nrn_pop1';
    model_prefix = sprintf('%s_%s_%gms_ls_%g_E_%s_P_%s',model_prefix_pre,nrn_model_ver,dur,...
                            layer_set_num,Efield_name,nrn_pop);
end
% Optional settings
in.data_fold = fullfile(mat_dir,'nrn_sim_data');
in.cell_ids = {1:5;6:10;11:15;16:20;21:25};
in.shift_dir = [-50,0,0]; 
in.plot_vertices = 1; % set to 1 plot patch FaceColors with 'interp'
% Data processing settings
in.comp_func = 'max'; % apply within cell (across compartments), e.g. 'max' or 'maxabs'
in.comp_type = 'soma'; % e.g. 'axon' or 'dendrite'
in.sectype = 'all'; % e.g. 'all' or 'terminal'
in.cell_func = 'median'; % function to apply across cells
in.norm_mode = 'none';
% Plot settings
in.ax = []; % axis handle
in.fig_units = 'centimeters';
in.fig_size = [50.8 28.36];
in.save_fig = 0;
in.fig_name = sprintf('pol_%s',model_prefix);
in.fig_fold = fullfile(mat_dir,'figures');
in.fig_format = '-png';
in.ax_view = [-89.2 45]; % [-89.2 70.8]
in.z_lims = [60 90.85]; 
in.lt_pos = [-411 -807 836]; % [170.3353 52.0246 1.2195e3]
in.clims = [];
in.cmap = 'bwr';
in.cbar_on = 0; 
% Input data rather than loading
in.deltaVms = []; 
in.cell_model_names = []; 
in = sl.in.processVarargin(in,varargin);
%% Load data
layers = loadLayers(layer_set_num);
if isempty(in.deltaVms)
    data_struct = load(fullfile(in.data_fold,[model_prefix '.mat']));
    deltaVms = data_struct.deltaVms;
    cell_model_names = data_struct.cell_model_names;
else
    deltaVms = in.deltaVms;
    cell_model_names = in.cell_model_names; 
    if ischar(cell_model_names)
       cell_model_names = {cell_model_names};  
    end
end
%% Extract data to plot
plot_cell_ids = [in.cell_ids{:}];
[~,~,data_inds] = intersect(cellModelNames(plot_cell_ids),cell_model_names);
deltaVms = deltaVms(data_inds);
cell_model_names = cell_model_names(data_inds);
%% Filter polarization values
if ~isempty(in.comp_func) && ~strcmp(in.comp_func,'none')
    deltaVms = processNrnCompData(deltaVms,plot_cell_ids,...
                            nrn_model_ver,in.comp_func,in.comp_type,in.sectype);                        
end
%% Plot
opts.func = in.cell_func;
opts.norm_mode = in.norm_mode;
opts.shift_dir = in.shift_dir;
opts.plot_vertices = in.plot_vertices;
opts.ax = in.ax;
pol_data = plotDataLayers(layers,deltaVms,cell_model_names,...
            in.cell_ids,opts);
if isempty(in.ax)
    ax = gca;
else
    ax = in.ax;
end
fig = ax.Parent;
fig.Units = in.fig_units;
fig.Position(3:4) = in.fig_size;
if ~isempty(in.clims)
    caxis(ax,in.clims);
end
if ischar(in.cmap)
    if strcmp(in.cmap,'bwr')
        colormap(ax,bluewhitered(1000));
    elseif strcmp(in.cmap,'rwb')
        colormap(ax,redwhiteblue(1000));
    end
else
    colormap(in.cmap);
end
view(ax,in.ax_view);
% change light
% ax = gca;
ax.Children(1).Position = in.lt_pos;
ax.Children(1).Style = 'local';
% cut off sulcus
if ~isempty(in.z_lims)
   ax.ZLim = in.z_lims;
end
% add colorbar
if in.cbar_on
    colorbar(ax,'FontSize',16);
end
%% Save figs
if in.save_fig
    printFig(fig,in.fig_fold,in.fig_name,'formats',{'fig','png'},'resolutions',{[],'-r250'});
    % savefig(fig,fullfile(in.fig_fold,[in.fig_name '.fig']));
    % export_fig(fig,fullfile(in.fig_fold,in.fig_name),in.fig_format,'-cmyk','-r250');
    % fprintf('Saved fig files to %s\n',fullfile(in.fig_fold,in.fig_name));
end
end
