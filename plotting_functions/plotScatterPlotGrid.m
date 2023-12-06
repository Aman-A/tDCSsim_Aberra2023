function [rsq_all,p_all,med_ape_all,mane2_all,mape_all] = plotScatterPlotGrid(x_data,...
                                                            y_data,varargin)
%PLOTSCATTERPLOTGRID ... 
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
in.save_fig = 0; 
in.fig_fold = '.';
in.fig_name = ''; 
in.fig_units = 'inches';
in.fig_size = [5 4]; 
in.font_size = 8;
in.font_name = 'Arial';
in.msize = 1; 
in.mcol = 0.6*ones(1,3);
in.fit_col = 'r';
in.unity_col = 'k';
in.y_lim = []; 
in.x_lim = []; 
in.add_xlabel = 1; 
in.ax_lw = 0.5; 
in.lw = 1; % for unity line
in.alpha = 0.005; 
in.median_in_layer = 0; % median across clones (1) or all points
in = sl.in.processVarargin(in,varargin); 
%% Plot
assert(length(x_data) == length(y_data),'Must have equal number of matrices in x_data and y_data'); 
num_layers = length(x_data); 
fig = figure('Color','w'); 
fig.Units = in.fig_units; 
fig.Position(3:4) = in.fig_size; 
rsq_all = zeros(num_layers,1); 
p_all = zeros(num_layers,1);
med_ape_all = zeros(num_layers,1); 
mane2_all = zeros(num_layers,1); 
mape_all = zeros(num_layers,1); 
all_sig = 1; 
for i=1:num_layers
    xi = x_data{i};
    yi = y_data{i};   
    if in.median_in_layer
        xi = nanmedian(xi,2);
        yi = nanmedian(yi,2); 
    end
    % 
    axi = subplot_tight(num_layers,1,i); 
%     axi = subplot(num_layers,1,i); 
    if isempty(xi)
        axis(axi,'off'); 
    else
        % Fits
        nan_inds = isnan(xi) | isnan(yi); % exclude nans from regressions to avoid nan output
        xi_fit = xi(~nan_inds); yi_fit = yi(~nan_inds); 
        statsi = calc_errs(xi_fit,yi_fit);
        med_ape = statsi.medape;        
        mane2 = statsi.mane2;
        lin_fit = polyfit(xi_fit,yi_fit,1);         
%         [R,p] = corrcoef(xi_fit,yi_fit); 
        rsq_all(i) = statsi.rsquared;
        p_all(i) = statsi.p_value;
        med_ape_all(i) = med_ape; 
        mane2_all(i) = mane2; 
        mape_all(i) = statsi.mape; 
%         rsq_all(i) = R(1,2)^2; 
%         p_all(i) = p(1,2); 
        if p_all(i) > in.alpha
            all_sig = 0; 
           fprintf('  L%g not significant, p = %f\n',i,p_all(i));  
        end   
        
        
%         p_err = 100*(xi_fit-yi_fit)./xi_fit; %percent error
%         med_ape = median(abs(p_err)); % median abs percent error
        fprintf('  L%g R^2 = %.3f, MedAPE: %.3f %%, mane2: %.3f\n',i,rsq_all(i),med_ape,mane2); 
        % Plot
        sc_lines = plot(axi,xi,yi,'.','MarkerSize',in.msize); hold on;    
        if length(in.mcol) == 1 || (isrow(in.mcol) && length(in.mcol) == 3)
            [sc_lines.Color] = deal(in.mcol);
        elseif ismatrix(in.mcol) && size(in.mcol,2) == 3
            for jj = 1:length(sc_lines)
    %            axi.Children(length(axi.Children)+1-jj).Color = in.mcol(jj,:);  
                sc_lines(jj).Color = in.mcol(jj,:); 
            end
        end
    %     if i == num_layers && in.add_xlabel
    %         xlabel('FEM E'); 
    %     end
    %     ylabel('Uniform E');
        axi.FontSize = in.font_size; 
        axi.FontName = in.font_name;         
        plot(axi,[min(xi,[],'all'),max(xi,[],'all')],[min(xi,[],'all'),max(xi,[],'all')],...
             '--','Color',in.unity_col,'LineWidth',in.lw); % add unity line
        plot(axi,[min(xi,[],'all'),max(xi,[],'all')],...
            polyval(lin_fit,[min(xi,[],'all'),max(xi,[],'all')]),...
            'Color',in.fit_col,'LineWidth',in.lw);    
        if isempty(in.x_lim)
            axi.XLim = [min(xi,[],'all'),max(xi,[],'all')];
        else
           axi.XLim = in.x_lim; 
        end
    %     axis(axi,'tight','square'); 
        axis(axi,'square'); 
        if isempty(in.y_lim)
            axi.YLim = axi.XLim; 
        else
            axi.YLim = in.y_lim;         
        end
        axi.LineWidth = in.ax_lw; 
        axi.YColor = 'k';
        axi.XColor = 'k';
        axi.Box = 'off';    
    %     axi.XTick = axi.YTick; 
        if min(xi(:),[],'all') < 0
            x_ticks = [min(xi(:),[],'all')*0.75,0,max(xi(:),[],'all')*0.75];
        else
            x_ticks = axi.XTick; 
        end
    %     x_ticks = [quantile(xi(:),0.02),0,quantile(xi(:),0.98)];
        if abs(max(xi(:))) < 0.01
            x_label = strsplit(num2str(x_ticks,'%.3f\n'));
            if any(x_ticks == 0)
                x_label{x_ticks == 0} = '0';
            end
        else
            x_label = strsplit(num2str(x_ticks,'%.2f\n'));
            if any(x_ticks == 0)
                x_label{x_ticks == 0} = '0';
            end
        end
        axi.XTick = x_ticks; 
        axi.XTickLabel = x_label;
        axi.YTick = x_ticks;
        axi.YTickLabel = x_label; 
    end
end
if all_sig
   fprintf('All regressions significant p < %f\n',in.alpha);  
end
if in.save_fig && ~isempty(in.fig_name)
    fig_name = in.fig_name;
    if in.median_in_layer
       fig_name = [fig_name '_med'];  
    end
    printFig(fig,in.fig_fold,fig_name,'formats',{'fig','png'},...
                'resolutions',{'','-r600'});    
end
end