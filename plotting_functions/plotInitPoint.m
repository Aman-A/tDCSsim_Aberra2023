% Enter cell_id, array of thetas and phis from full threshold simulation,
% and cell array of initial spike secnames corresponding to the
% orientations in thetas, phis. Plots AP initiation point on cell
% morphology using compartment coordinates

function plotInitPoint(ax1,cell_id,data,thetas,phis,init_inds,...
         arrow_size,expand_ax,xlimits,zlimits,mark_size,nrn_model_ver,varargin)       
    in.plot_ind = 'min';
    in.comp_data = [];
    in.lw = 0.5; % line width for plotCellLines
    in.plot_comp_type = 'all';
    in.cell_models_file = 'cell_models';
    in = sl.in.processVarargin(in,varargin);
    if isempty(init_inds)
        plot_sec = 0;
    else
        plot_sec = 1;
    end
    if ischar(in.plot_ind)
        if strcmp(in.plot_ind,'min')
            [~,plot_ind] = min(data);
        elseif strcmp(in.plot_ind,'max')        
            [~,plot_ind] = max(data);    
        elseif strcmp(in.plot_ind,'maxabs')        
            [~,plot_ind] = max(abs(data));    
        end
    elseif isnumeric(in.plot_ind)
        plot_ind = in.plot_ind;
    end
    theta_plot = thetas(plot_ind); phi_plot = phis(plot_ind);  
%     fprintf('Plotting theta %d phi %d\n',theta_plot,phi_plot); 
    % Get cell coordinates        
    cell_data0 = loadCellData(cell_id,nrn_model_ver,...
                        {'C','comp_types','parent_inds','sectypes','diams'},...
                        'cell_models_file',in.cell_models_file); % original cell
    if ~strcmp(in.plot_comp_type,'all')
       cell_data = extractCellSubMorph(cell_data0,in.plot_comp_type); 
       if size(in.comp_data,2) ~= size(cell_data.C,1) % probably recorded all comps 
            comp_inds = find(getCellCompInds(cell_data0,in.plot_comp_type,'all'));            
       else
           comp_inds = find(true(size(cell_data.C,1))); 
       end
    else
        cell_data = cell_data0; 
        comp_inds = find(true(size(cell_data0.C,1),1)); 
    end
    C = cell_data.C; 
    % Plot     
    [Ex,Ey,Ez] = angle2vec(theta_plot,phi_plot,1);
    axes(ax1); % set ax1 as gca
    if isempty(in.comp_data)
        plotCellLines('cell_data',cell_data,'lw',in.lw);
    else
        plotCellLines('cell_data',cell_data,'vals',in.comp_data(plot_ind,comp_inds)',...
                      'lw',in.lw,'plot_soma_sphere',0);
        colormap(ax1,bluewhitered(1000)); % polarization 
    end
    view(ax1,[0 -1 0]);                
    ax1.DataAspectRatio = [1 1 1];    
    clims = caxis(ax1);      
    caxis(ax1,'manual'); % fix colormap scaling     
    hold(ax1,'on');
    grid_mult = 1;
    [X,Y,~] = meshgrid(grid_mult*linspace(min(C(:,1)),max(C(:,1)),3),grid_mult*linspace(min(C(:,2)),max(C(:,2)),3),grid_mult*linspace(min(C(:,3)),max(C(:,3)),3));
    %[X,Y,Z] = meshgrid(grid_mult*linspace(-1365,1169,3),grid_mult*linspace(min(C(:,2)),max(C(:,2)),3),grid_mult*linspace(-1170,390.2,3));
%     U = Ex*ones(size(X)); V = Ey*ones(size(X)); W = Ez*ones(size(X));
     
    if plot_sec
        sec_ind = init_inds(plot_ind); 
        % Plot AP initiation point
        if isempty(sec_ind)
           disp('Min threshold not found'); 
        end        
        plot3(ax1,C(sec_ind,1),min(C(:,2)),C(sec_ind,3),'Marker','p','MarkerSize',mark_size,'MarkerFaceColor','white','MarkerEdgeColor','black',...
            'LineWidth',0.5);
    end
    sm = 1;

    if ~strcmp(xlimits,'auto')
        ax1.XLim = xlimits; 
    end
    if ~strcmp(zlimits,'auto')
        ax1.ZLim = zlimits;
    end    
%     [~,max_dim] = max([Ex,Ey,Ez]); 
%     arrow_size = 0.6*range(C(:,max_dim));
    if Ex < 0
        x_shift = -Ex*arrow_size;
    else
       x_shift = 0; 
    end
    if Ez > 0
        z_shift = -Ez*arrow_size;
    else
        z_shift = 0;
    end
    quiver3D([ax1.XLim(1)*0.8+x_shift,0,0.8*ax1.ZLim(2)+z_shift],arrow_size*[Ex Ey Ez]); 
    camlight head; 
    % Plot scale bar in bottom left of figure
    plot3(ax1,[sm*ax1.XLim(1),sm*ax1.XLim(1),sm*ax1.XLim(1)+250],[min(min(min(Y))),min(min(min(Y))),min(min(min(Y)))],[ax1.ZLim(1)+250,ax1.ZLim(1),ax1.ZLim(1)],...
        'Color','k','LineWidth',3);
    hold off;
    % reset color limits to original limits
    ax1.CLim = clims;    
    if expand_ax   
        ax1.ZLim = [ax1.ZLim(1)*1.1 ax1.ZLim(2)]; 
    end 
    
end