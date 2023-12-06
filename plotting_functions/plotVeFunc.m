% input an axis handle, the model_prefix, cell_model_name, vector of
% thresholdEs, vector of theta and phi angles, cell array of section at
% which spike initiated, and corresponding vector of section types
% plots the cell morphology with a vector indicating direction of E-field
% for minimum threshold and a marker on the point of initiation, with the
% marker type corresponding to the section type
%
% Aman Aberra

function plotVeFunc(ax1,cellid,data,thetas,phis,init_inds,nrn_model_ver,varargin)    
    in.plot_ind = 'min';  
    in.comp_data = []; 
    in.lw = 0.5;
    in.plot_comp_type = 'all';
    in.cell_models_file = 'cell_models';
    in = sl.in.processVarargin(in,varargin); 
    expand_ax = 0;
%     x_arrow_shift = -200;
%     z_arrow_shift = 170;    
    cell_model_name = cellModelNames(cellid,...
                                    'cell_models_file',in.cell_models_file);
    [~,~,layer_mtype,~,clone_num] = cellNameParser(cell_model_name);
    mout = outputMorphParams(nrn_model_ver,struct('cell_id',cellid,'cell_model_name',cell_model_name),...
                    'print_output',0,'cell_models_file','cell_models_all');
    if mout.replace_axon == 1 || strcmp(in.plot_comp_type,'soma_dendrite')
       includes_axon = 0; % use tighter axis limits 
    else
       includes_axon = 1;
    end
    if strcmp(layer_mtype,'L1_NGC-DA')
        if includes_axon
            xlimits = [-400 400]; zlimits = [-200 200];
            arrow_size = 316;
        else
            xlimits = [-250 200]; zlimits = [-200 100];
            arrow_size = 200;
        end
    elseif strcmp(layer_mtype,'L23_PC')    
        if includes_axon
            xlimits = [-1500 1100]; 
            if clone_num == 2
                zlimits = [-800 750]; 
            else
%                 zlimits = [-1200 350]; 
                zlimits = [-1200 500]; 
            end
            arrow_size = 1000;
        else
            xlimits = [-400 400]; 
            if clone_num == 2
                zlimits = [-400 600]; 
            else
                zlimits = [-400 600]; 
            end
            arrow_size = 300;
        end
        
    elseif strcmp(layer_mtype,'L4_LBC')            
        if includes_axon
            xlimits = [-700 500]; zlimits = [-800 800]; 
            arrow_size = 620; 
        else
            xlimits = [-400 400]; zlimits = [-600 300]; 
            arrow_size = 300; 
        end
    elseif strcmp(layer_mtype,'L5_TTPC2')
        if includes_axon
            xlimits = [-1000 1000]; zlimits = [-1000 1500]; 
            arrow_size = 850;
        else
            xlimits = [-400 400]; zlimits = [-400 1500]; 
            arrow_size = 500;
        end
    elseif strcmp(layer_mtype,'L6_TPC_L4') || strcmp(layer_mtype,'L6_TPC_L1')
        if includes_axon
            xlimits = [-1000 1000]; 
            if clone_num == 1
                zlimits = [-720 1500];
            else
                zlimits = [-750 1500];
            end
            arrow_size=600;
            expand_ax = 1;
%             xlimits = [-1000 1000]; zlimits = [-1000 1500];
%             arrow_size = 850;
        else
            xlimits = [-300 500]; 
            if clone_num == 1
                zlimits = [-400 1500];
            else
                zlimits = [-400 1500];
            end
            arrow_size=400;
            expand_ax = 1;
        end
    else        
        xlimits = [-1000 1000]; zlimits = [-1000 1500]; 
        arrow_size = 850;        
    end           
    mark_size = 10;    
    plotInitPoint(ax1,cellid,data,thetas,phis,init_inds,arrow_size,...
                expand_ax,xlimits,zlimits,mark_size,nrn_model_ver,in)
end