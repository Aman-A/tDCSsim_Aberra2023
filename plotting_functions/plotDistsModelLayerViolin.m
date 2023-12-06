function data_roi = plotDistsModelLayerViolin(data_all,ROIi,cell_ids,...
                                              cell_model_names,varargin)
%PLOTDISTSMODELLAYER Plots boxplots for models across model, subgrouped by layer
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
in.colors = jet(length(cell_ids)); 
in.violin_lw = 2;
in.violin_face_alpha = 0.5; 
in.plot_pts = 0; 
in.pts_col = 0.8*ones(1,3); 
in.pt_markerSize = 1; 
in = sl.in.processVarargin(in,varargin); 
num_models = length(data_all); 
num_layers = length(cell_ids); 
%%
% combine data within each layer
data = cell(num_layers,num_models); % rows - layer, columns - model
for j = 1:num_models
    for i = 1:num_layers
        if ~isempty(cell_ids{i})
            cell_model_names_i = cellModelNames(cell_ids{i}); % cell names in layer
            [~,~,data_inds] = intersect(cell_model_names_i,cell_model_names); % get indices of layer cells in threshEs
            if ~isempty(data_inds)
                data_layer =  cell2mat(data_all{j}(data_inds));
                data{i,j} = data_layer; % combine all thresholds in layer
            else
                data{i,j} = 1e6; % placeholder large number to allow UniformOutput in cellfun below, no cells in this layer
            end
        else
           data{i,j} = 1e6;  
        end
    end
end
% Extract data within desired region
layer_names = strcat('L',strsplit(num2str(1:num_layers)))';
layer_ids_all = {}; 
data_roi = data; 
data_plot = cell(size(data_roi));
for j = 1:num_models 
    for i = 1:num_layers  % layer specific ROI indices              
        if size(data_roi{i,j},1) == size(ROIi{i},1)
            data_roi{i,j} = data{i,j}(ROIi{i},:);
            data_plot{i,j} = data_roi{i,j}(:); % turn into 1D array
            layer_ids_all = [layer_ids_all; repmat(strcat(layer_names(i), '_m',num2str(j)),length(data_plot{i}),1)];                   
        end
    end
end
% threshEs_all = threshEs_all'; 
% data_plot = cell2mat(data_plot(:)); % collapse into single vector
%% Plot
% boxplot(data_plot,layer_ids_all,'Symbol','','Width',0.5); 
% boxplot(data_plot,layer_ids_all,'Width',0.5,'Whisker',whisker); 
x_vals = zeros(size(data_plot,1),size(data_plot,2)); 
for j = 1:size(data_plot,2)
    for i = 1:size(data_plot,1)
        x_vals(i,j) = i + size(data_plot,1)*(j-1);
        if in.plot_pts && ~isempty(data_plot{i,j})
%             plot(x_vals(i,j)*ones(size(data_plot{i,j})),data_plot{i,j},'k.','MarkerSize',4)
            plotSpread(data_plot{i,j},...
                       'xValues',x_vals(i,j),...
                       'distributionColors',0.8*ones(1,3),...
                       'markerSize',in.pt_markerSize); hold on;
        end
        violin_func(data_plot{i,j},'x',x_vals(i,j),'facecolor',in.colors(i,:),...
                'mc',[],'plotlegend',0,'lw',in.violin_lw,'facealpha',in.violin_face_alpha);   
        hold on;           
    end
end
ax = gca;
ax.XLim = [min(x_vals,[],'all')-1,max(x_vals,[],'all')+1];
yaxlim = ax.YLim*10;
num_plot_layers = num_layers; 
% num_plot_layers = sum(cellfun(@length,data_roi(:,1))>1);
if num_models > 1 % plot dividing lines    
    ax.XTick = ((num_plot_layers+1)/2:num_plot_layers:(num_models*num_plot_layers-1)+1);    
%     x_ticks = x_vals(plot_layers,:);
%     ax.XTick = x_ticks(:); 
    ax.XTickLabel = strcat('M',strsplit(num2str(1:num_models)))';
    if strcmp(ax.YScale,'log')        
        y_min = 1e-3*min(data_plot);
    else        
        y_min = yaxlim(1); 
        if y_min < 0
           y_min = y_min*3; 
        else
           y_min = y_min*0.4;
        end
    end
    plot(repmat(num_plot_layers:num_plot_layers:num_plot_layers*(num_models-1),2,1)+0.5,...
                [y_min*ones(1,num_models-1);3*yaxlim(2)*ones(1,num_models-1)],...
                'Color','k','LineWidth',ax.LineWidth)
%     leg_objs = bp.Children(num_models*num_plot_layers+1:2*num_models*num_plot_layers);
    leg_objs = ax.Children(2:num_models*num_plot_layers); 
    leg_objs = leg_objs(num_plot_layers:-1:1); 
    legend(leg_objs,layer_names(~cellfun(@isempty,cell_ids)));
end
end
