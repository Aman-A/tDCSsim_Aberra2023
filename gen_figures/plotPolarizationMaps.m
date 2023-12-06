%% Plot grid of somatic polarization-direction maps for all model neurons
% Set fig_name to get predefined settings replicating figure panels
% or set to 'cust' to define settings below

fig_name = 'FigS5A'; % valid options: 'Fig2D', 'FigS2','FigS3', 
                    % 'FigS4', or 'FigS5A'
save_fig = 0; % set to 1 to save
if strcmp(fig_name,'cust')
    nrn_model_ver = 'maxH'; % 'maxH' for myelinated, 'umaxH' for unmyelinated
    comp_func = 'max'; % function to apply within cell (across compartments)
                       % 'max', 'maxabs','min', etc.
    comp_type = 'soma'; % 'soma', 'axon', or 'dendrite'
    sectype = 'all'; % 'all' or 'terminal'
    plot_comp_type = 'all'; % 'all'  to visualize whole morphology or 
                            % 'soma_dendrite' or 'soma_axon' to visualize soma
                            % + dendrite or soma + axon only, respectively
else
    switch fig_name
        case 'Fig2D' % somatic polarization, myelinated axons
            nrn_model_ver = 'maxH';
            comp_func = 'max'; 
            comp_type = 'soma';
            sectype = 'all'; 
            plot_comp_type = 'all';
            reorder = 1;
        case 'FigS2' % peak axon terminal polarization, myelinated axons
            nrn_model_ver = 'maxH';
            comp_func = 'maxabs'; 
            comp_type = 'axon'; 
            sectype = 'terminal';
            plot_comp_type = 'soma_axon';
            reorder = 1;
        case 'FigS3' % peak basal dendritic terminal polarization, myelinated axons
            nrn_model_ver = 'maxH';
            comp_func = 'maxabs'; 
            comp_type = 'basal'; 
            sectype = 'terminal';
            plot_comp_type = 'soma_axon';
            reorder = 1;
        case 'FigS4' % peak apical dendritic terminal polarization, myelinated axons
            nrn_model_ver = 'maxH';
            comp_func = 'maxabs';
            comp_type = 'apical';
            sectype = 'terminal';
            plot_comp_type = 'soma_axon';
            reorder = 1;
        case 'FigS5A' % peak somatic terminal polarization, unmyelinated axons
            nrn_model_ver = 'umaxH';
            comp_func = 'max'; 
            comp_type = 'soma';
            sectype = 'all'; 
            plot_comp_type = 'all';
            % use same ordering as Fig 2D (myelinated)
            reorder = [2 8 12 20 23 4 10 13 16 24 5 9 11 17 25 3 7 15 19 22 1 6 14 18 21];
    end
end
plotPol_MWgrid_morphs(save_fig,nrn_model_ver,comp_func,comp_type,sectype,...
                      plot_comp_type,'order_by_maxpol',reorder)                        


