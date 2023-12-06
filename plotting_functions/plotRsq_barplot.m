function plotRsq_barplot(data_all1,data_all2,data_all3,...
                         sig_all,save_fig,varargin)
% can enter empty inputs for data_all2/data_all3, same size as data_all1
% sig_all: num_layers x num_comp_types x num_data                      
in.fig_size = [6.3 2.2]; 
in.font_name = 'Arial';
in.font_size = 10; 
in.cell_labels = {'L1','L2/3','L4','L5','L6'}; 
in.title_on = 0; 
in.bar_cols = [27,158,119;217,95,2;117,112,179]*(1/255); 
in.x_shifts = [-2/9,0,2/9];
in.bar_width = 0.2;
in.comp_types = {'soma','axon','basal','apical'}; 
in.fig_name = 'Rsq_bars';
in.fig_fold = './';
in.YTick = 0:0.25:1; % default for R2
in.y_lim = []; 
in.y_lim_same = 0; % only necessary if y_lim not being set explicitly
in = sl.in.processVarargin(in,varargin); 
%%
num_comp_types = length(data_all1); 
num_layers = size(data_all1{1},1);
%% Plot
fig = figure('Color','w');
fig.Units = 'inches';
fig.Position(3:4) = in.fig_size; 
y_lims_all = zeros(num_comp_types,2); 
for c = 1:num_comp_types    
    axc = subplot_tight(1,num_comp_types,c,[0.15 0.03]); % [ vert horz] margins
    all_axes{c} = axc; 
    datai = [data_all1{c}]; % '|E| E_n, E_t
    if ~isempty(data_all2{c})
       datai = [datai,data_all2{c}]; 
    end
    if ~isempty(data_all3{c})
       datai = [datai,data_all3{c}]; 
    end
    for i = 1:num_layers        
        if ~isempty(in.bar_cols)
           for j = 1:size(datai,2)   
               if size(datai,2) == 1 && size(in.bar_cols,1) == num_layers
                  bar_colij = in.bar_cols(i,:); % color by layer
               else
                   bar_colij = in.bar_cols(j,:); % color by data within layer
               end
               % [bottom left, top left, top right, bottom right] (CW)
               xcoords = [i+in.x_shifts(j) - in.bar_width/2,i+in.x_shifts(j) - in.bar_width/2,...
                          i+in.x_shifts(j) + in.bar_width/2,i+in.x_shifts(j) + in.bar_width/2];
               ycoords = [0 datai(i,j) datai(i,j) 0]; 
               b = patch(xcoords,ycoords,...
                         bar_colij,'EdgeColor','none');
               if isequal(bar_colij,[1 1 1])
                    b.EdgeColor = in.bar_cols(j-1,:);
               end
%               b = bar(i+in.x_shifts(j),datai(i,j),'BarWidth',in.bar_width); hold on;
%               b.FaceColor = in.bar_cols(j,:);  
%               b.EdgeColor = 'none';
               if sig_all(i,c,j) == 0 
                  hatchfill2(b,'Single','HatchColor','k');  
               end
           end
        end
    end
    box(axc,'off'); 
    axc.YGrid = 'on';
    axc.FontName = in.font_name;
    axc.FontSize = in.font_size; 
    axc.XTick = 1:num_layers; 
    axc.XTickLabel = in.cell_labels;
    axc.XColor = 'k';
    axc.YColor = 'k';
    axc.XLim = [0.5 num_layers+0.5]; 
    if ~isempty(in.y_lim)
       axc.YLim = in.y_lim;  
    end
%     if i == 1
%        ylabel('R^{2}') 
%     end    
    if in.title_on
       title(in.comp_types{c});  
    end
    if ~isempty(in.YTick)
        axc.YTick = in.YTick;        
        if c > 1
            axc.YTickLabel = {};
        end
    end    
    axc.Position(1) = axc.Position(1)+0.025;
    y_lims_all(c,:) = axc.YLim; 
end
if in.y_lim_same
    y_lim = [min(y_lims_all(:,1)),max(y_lims_all(:,2))]; 
   for c = 1:num_comp_types 
       all_axes{c}.YLim = y_lim; % set all to min/max across panels       
       if c == 1
           yticks1 = all_axes{c}.YTick; 
       else           
           all_axes{c}.YTick = yticks1; 
           all_axes{c}.YTickLabel = {}; % can remove since all have same limits
        end
   end
end
if save_fig    
    printFig(fig,in.fig_fold,in.fig_name,'formats',{'fig','eps'});
end
end
                     
                     