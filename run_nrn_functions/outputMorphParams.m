% input nrn_model_ver and current params_struct m, or no params struct
% make sure min has cell_id field for nrn_model_ver's that have cell_id
% specific settings (e.g. 'maxH')
function mout = outputMorphParams(nrn_model_ver,min,varargin)
in.print_output = 1; 
in.cell_models_file = 'cell_models';
in = sl.in.processVarargin(in,varargin); 
existing_nrn_model_vers = {'umaxr','max4rP60c','maxH','maxmH','umaxH','maxM','maxM_giganto','maxH_rax1',...
                            'max4rP60c_rax1','umaxr_rax1','maxH2','maxH3',...
                            'max4rP60c2','max4rP60c3'}; 
nrn_model_vers_list_str = sprintf('%s\n',existing_nrn_model_vers{:}); 
if nargin == 0
   fprintf('Enter one of the following nrn_model_ver names:\n');
   fprintf(nrn_model_vers_list_str); 
   mout = []; 
   return
end
if nargin == 1
   min = struct; % initialize to empty struct
end
mout = min;
mout.nrn_model_ver = nrn_model_ver;    
if ~isfield(min,'cell_models_file')
   mout.cell_models_file = in.cell_models_file;     
end
err_str = [sprintf('Input %s not recognized, existing nrn_model_vers:\n',nrn_model_ver),...
            nrn_model_vers_list_str]; 
assert(any(strcmp(nrn_model_ver,existing_nrn_model_vers)),...
        sprintf(err_str)); 
    
if strcmp(nrn_model_ver,'umaxr')
    % P14 rat, unmyelinated, no scaling
    mout.myelinate_axon = 0;
    mout.prune_axon = 0;
    mout.replace_axon = 0;
    mout.scale_axon_diam = 1;
    mout.scale_apic_diam = 1;
    mout.scale_basal_diam = 1;
    mout.scale_soma_area = 1;
    mout.scale_basal_L = 1;
    mout.min_myelinD = 0.2; 
    mout.min_myelinL = 20; 
    mout.max_myelin_order = 0;   
elseif strcmp(nrn_model_ver,'max4rP60c')
    % adult, rat with myelinated axon, using Zhu 2000 scaling (Aberra 2018)
    mout.myelinate_axon = 1;
    mout.prune_axon = 0;
    if isL5_6PC(mout)
        mout.replace_axon = 2; % disable main axon terminal
    else
        mout.replace_axon = 0;
    end
    mout.scale_axon_diam = 1.322; % using Zhu 2000 L5 pc soma scaling
    mout.scale_apic_diam = 1.248; % using Romand 2011 L5 PC growth
    mout.scale_basal_diam = 1.133; % using Romand 2011 L5 PC growth
    mout.scale_soma_area = 1.322; % using Zhu 2000 L5 pc soma scaling
    mout.scale_basal_L = 1.17; % using Romand 2011 L5 PC growth
    mout.min_myelinD = 0.2; 
    mout.min_myelinL = 20; 
    mout.max_myelin_order = 0;
    mout.max_myelin_order_mode = 'relative'; 
    mout.INL_ratio = 100;
    mout.INL_ratio_term = 70;
    mout.nodeL = 1;
    mout.min_PMAS = 50; 
    mout.myelinL_error = 0.1; 
    mout.nodeL_error = 0.1; 
    mout.min_pretermL = 0; 
    mout.myelin_cm_mode = 'constant'; % sets to 0.02 uF/cm2
    mout.myelin_gpas_mode = 'constant'; % sets to 1/1.125e6 
elseif strcmp(nrn_model_ver,'maxH') || strcmp(nrn_model_ver,'maxmH')
    % adult, human with myelinated axon using L2/3 PC rat BB: human Eyal
    % 2018 (Aberra 2018/Aberra 2020)
    mout.myelinate_axon = 1;
    mout.prune_axon = 0;  
    if isL5_6PC(mout)
        mout.replace_axon = 2; % disable main axon terminal
    else
        mout.replace_axon = 0;
    end
    mout.scale_axon_diam = 2.453; % use soma scaling
    mout.scale_apic_diam = 1.876;
    mout.scale_basal_diam = 1.946;
    mout.scale_soma_area = 2.453;
    mout.scale_basal_L = 1.17; % using Romand 2011 L5 PC growth
    mout.min_myelinD = 0.2;
    mout.min_myelinL = 20; 
    if strcmp(nrn_model_ver,'maxmH') % myelinate main axon  only
        mout.max_myelin_order = 1;
        mout.max_myelin_order_mode = 'absolute';  
    else
        mout.max_myelin_order = 0;
        mout.max_myelin_order_mode = 'relative'; 
    end
    mout.INL_ratio = 100;
    mout.INL_ratio_term = 70;
    mout.nodeL = 1;
    mout.min_PMAS = 50; 
    mout.myelinL_error = 0.1; 
    mout.nodeL_error = 0.1; 
    mout.min_pretermL = 0; 
    mout.myelin_cm_mode = 'constant'; % sets to 0.02 uF/cm2
    mout.myelin_gpas_mode = 'constant'; % sets to 1/1.125e6 
elseif strcmp(nrn_model_ver,'maxM')
    % monkey with myelinated axon using 2 different sets of scaling factors
    mout.myelinate_axon = 1;
    mout.prune_axon = 0;  
    if isL5_6PC(mout)
        mout.replace_axon = 2; % disable main axon terminal
        mout.scale_axon_diam = 1.375; % use soma scaling
        mout.scale_apic_diam = 0.987;
        mout.scale_basal_diam = 1;
        mout.scale_soma_area = 1.14;
    else
        mout.replace_axon = 0;
        mout.scale_axon_diam = 1.59; % use soma scaling
        mout.scale_apic_diam = 1.10;
        mout.scale_basal_diam = 1.1;
        mout.scale_soma_area = 1.59;
    end
    
    mout.scale_basal_L = 1.17; % using Romand 2011 L5 PC growth
    mout.min_myelinD = 0.2;
    mout.min_myelinL = 20; 
    mout.max_myelin_order = 0;
    mout.max_myelin_order_mode = 'relative'; 
    mout.INL_ratio = 100;
    mout.INL_ratio_term = 70;
    mout.nodeL = 1;
    mout.min_PMAS = 50; 
    mout.myelinL_error = 0.1; 
    mout.nodeL_error = 0.1; 
    mout.min_pretermL = 0; 
    mout.myelin_cm_mode = 'constant'; % sets to 0.02 uF/cm2
    mout.myelin_gpas_mode = 'constant'; % sets to 1/1.125e6 

elseif strcmp(nrn_model_ver,'maxM_giganto')
    % monkey with myelinated axon using gigantopyramidal cells scaling
    % factors
    mout.myelinate_axon = 1;
    mout.prune_axon = 0;  
    if isL5_6PC(mout)
        mout.replace_axon = 2; % disable main axon terminal
        mout.scale_axon_diam = 3.27; % use soma scaling
        mout.scale_apic_diam = 1.90;
        mout.scale_basal_diam = 2.16;
        mout.scale_soma_area = 2.07;
    else
        mout.replace_axon = 0;
        mout.scale_axon_diam = 1.59; % use soma scaling
        mout.scale_apic_diam = 1.10;
        mout.scale_basal_diam = 1.1;
        mout.scale_soma_area = 1.59;
    end
    
    mout.scale_basal_L = 1.17; % using Romand 2011 L5 PC growth
    mout.min_myelinD = 0.2;
    mout.min_myelinL = 20; 
    mout.max_myelin_order = 0;
    mout.max_myelin_order_mode = 'relative'; 
    mout.INL_ratio = 100;
    mout.INL_ratio_term = 70;
    mout.nodeL = 1;
    mout.min_PMAS = 50; 
    mout.myelinL_error = 0.1; 
    mout.nodeL_error = 0.1; 
    mout.min_pretermL = 0; 
    mout.myelin_cm_mode = 'constant'; % sets to 0.02 uF/cm2
    mout.myelin_gpas_mode = 'constant'; % sets to 1/1.125e6 
    
elseif strcmp(nrn_model_ver,'umaxH')
    % adult, human with unmyelinated axon using L2/3 PC rat BB: human Eyal 2018
    mout.myelinate_axon = 0;
    mout.prune_axon = 0;
    mout.replace_axon = 0;
    mout.scale_axon_diam = 2.453; 
    mout.scale_apic_diam = 1.876;
    mout.scale_basal_diam = 1.946;
    mout.scale_soma_area = 2.453;
    mout.scale_basal_L = 1.17;     
    
elseif strcmp(nrn_model_ver,'maxH_rax1')
    % scale soma/dendrite to adult, human, but replace axon with just AIS
    mout.myelinate_axon = 0;
    mout.prune_axon = 0;
    mout.replace_axon = 1;
    mout.scale_axon_diam = 2.453; 
    mout.scale_apic_diam = 1.876;
    mout.scale_basal_diam = 1.946;
    mout.scale_soma_area = 2.453;
    mout.scale_basal_L = 1.17; 
elseif strcmp(nrn_model_ver,'max4rP60c_rax1')
    % scale adult, rat with myelinated axon, using Zhu 2000, but replace
    % axon with just AIS
    mout.myelinate_axon = 0;
    mout.prune_axon = 0;
    mout.replace_axon = 1;
    mout.scale_axon_diam = 1.322; % using Zhu 2000 L5 pc soma scaling
    mout.scale_apic_diam = 1.248; % using Romand 2011 L5 PC growth
    mout.scale_basal_diam = 1.133; % using Romand 2011 L5 PC growth
    mout.scale_soma_area = 1.322; % using Zhu 2000 L5 pc soma scaling
    mout.scale_basal_L = 1.17; % using Romand 2011 L5 PC growth    
elseif strcmp(nrn_model_ver,'umaxr_rax1')
    % P14 rat, unmyelinated, no scaling, replace axon with just AIS - same
    % as original BB cell (temp was 34 degC in original models)
    mout.myelinate_axon = 0;
    mout.prune_axon = 0;
    mout.replace_axon = 0;
    mout.scale_axon_diam = 1;
    mout.scale_apic_diam = 1;
    mout.scale_basal_diam = 1;
    mout.scale_soma_area = 1;
    mout.scale_basal_L = 1;
    mout.replace_axon = 1; 
elseif strcmp(nrn_model_ver,'maxH2')
    % adult, human with myelinated axon using L2/3 PC rat BB: human Eyal 2018
    % version 2
    mout.myelinate_axon = 1;
    mout.prune_axon = 0;
    if isL5_6PC(mout)
        mout.replace_axon = 2; % disable main axon terminal
    else
        mout.replace_axon = 0;
    end
    mout.scale_axon_diam = 2.453; % use soma scaling
    mout.scale_apic_diam = 1.876;
    mout.scale_basal_diam = 1.946;
    mout.scale_soma_area = 2.453;
    mout.scale_basal_L = 1.17; % using Romand 2011 L5 PC growth
    mout.min_myelinD = 0.2;
    mout.min_myelinL = 20; % 15? (Stedehouder)
    mout.max_myelin_order = 10;
    mout.max_myelin_order_mode = 'absolute'; 
    mout.INL_ratio = 100;
    mout.INL_ratio_term = 70;
    mout.nodeL = 1; % 1.5 (Stedehouder)
    mout.min_PMAS = 50; 
    mout.myelinL_error = 0.1; 
    mout.nodeL_error = 0.1; 
    mout.min_pretermL = 30; 
    mout.myelin_cm_mode = 'variable'; % sets based on myelin thickness
    mout.myelin_gpas_mode = 'variable'; % sets based on myelin thickness
    mout.include_rpa = 0; % includes periaxonal resistivity
    mout.myelin_layers = 1; % 2 layers for myelin
elseif strcmp(nrn_model_ver,'maxH3')
     % adult, human with myelinated axon using L2/3 PC rat BB: human Eyal 2018
    % version 3
    mout.myelinate_axon = 1;
    mout.prune_axon = 0;
    if isL5_6PC(mout)
        mout.replace_axon = 2; % disable main axon terminal
    else
        mout.replace_axon = 0;
    end
    mout.scale_axon_diam = 2.453; % use soma scaling
    mout.scale_apic_diam = 1.876;
    mout.scale_basal_diam = 1.946;
    mout.scale_soma_area = 2.453;
    mout.scale_basal_L = 1.17; % using Romand 2011 L5 PC growth
    mout.min_myelinD = 0.2;
    if isPC(mout)
        mout.min_myelinL = 30; % Cohen et al. 2020 minimum internode observed
        mout.max_myelin_order = 10;
        mout.max_myelin_order_mode = 'absolute'; 
    else
        mout.min_myelinL = 20; % 15? (Stedehouder)
        mout.max_myelin_order = 10;
        mout.max_myelin_order_mode = 'absolute'; 
    end    
    mout.INL_ratio = 100;
    mout.INL_ratio_term = 70;
    mout.nodeL = 1; % 1.5 (Stedehouder)
    mout.min_PMAS = 50; 
    mout.myelinL_error = 0.1; 
    mout.nodeL_error = 0.1; 
    mout.min_pretermL = 30; 
    mout.myelin_cm_mode = 'variable'; % sets based on myelin thickness
    mout.myelin_gpas_mode = 'variable'; % sets based on myelin thickness
    mout.include_rpa = 0; % includes periaxonal resistivity
    mout.myelin_layers = 1; % 2 layers for myelin
elseif strcmp(nrn_model_ver,'max4rP60c2')
    % adult, rat with myelinated axon, using Zhu 2000 scaling
    % version 2
    mout.myelinate_axon = 1;
    mout.prune_axon = 0;
    if isL5_6PC(mout)
        mout.replace_axon = 2; % disable main axon terminal
    else
        mout.replace_axon = 0;
    end
    mout.scale_axon_diam = 1.322; % using Zhu 2000 L5 pc soma scaling
    mout.scale_apic_diam = 1.248; % using Romand 2011 L5 PC growth
    mout.scale_basal_diam = 1.133; % using Romand 2011 L5 PC growth
    mout.scale_soma_area = 1.322; % using Zhu 2000 L5 pc soma scaling
    mout.scale_basal_L = 1.17; % using Romand 2011 L5 PC growth
    mout.min_myelinD = 0.2; 
    mout.min_myelinL = 20; 
    mout.max_myelin_order = 10;
    mout.max_myelin_order_mode = 'absolute'; 
    mout.INL_ratio = 100;
    mout.INL_ratio_term = 70;
    mout.nodeL = 1;
    mout.min_PMAS = 50; 
    mout.myelinL_error = 0.1; 
    mout.nodeL_error = 0.1; 
    mout.min_pretermL = 30; 
    mout.myelin_cm_mode = 'variable'; % sets to 0.02 uF/cm2
    mout.myelin_gpas_mode = 'variable'; % sets to 0.02 uF/cm2
    mout.include_rpa = 0; % includes periaxonal resistivity
    mout.myelin_layers = 1; 
elseif strcmp(nrn_model_ver,'max4rP60c3')
    % adult, rat with myelinated axon, using Zhu 2000 scaling, Stehedouder
    % 2019 and Cohen 2020 parameters
    % version 3
    mout.myelinate_axon = 1;
    mout.prune_axon = 0;
    if isL5_6PC(mout)
        mout.replace_axon = 2; % disable main axon terminal
    else
        mout.replace_axon = 0;
    end
    mout.scale_axon_diam = 1.322; % using Zhu 2000 L5 pc soma scaling
    mout.scale_apic_diam = 1.248; % using Romand 2011 L5 PC growth
    mout.scale_basal_diam = 1.133; % using Romand 2011 L5 PC growth
    mout.scale_soma_area = 1.322; % using Zhu 2000 L5 pc soma scaling
    mout.scale_basal_L = 1.17; % using Romand 2011 L5 PC growth
    mout.min_myelinD = 0.33; % Stedehouder 2019     
    if isPC(mout)
        mout.min_myelinL = 30; % Cohen et al. 2020 minimum internode observed
    else
        mout.min_myelinL = 14; % Stedehouder 2019
    end
    mout.max_myelin_order = 10;
    mout.max_myelin_order_mode = 'absolute'; 
    mout.INL_ratio = 100;
    mout.INL_ratio_term = 70;
    mout.nodeL = 1;
    mout.min_PMAS = 50; 
    mout.myelinL_error = 0.1; 
    mout.nodeL_error = 0.1; 
    mout.min_pretermL = 30; 
    mout.myelin_cm_mode = 'variable'; % sets to 0.02 uF/cm2
    mout.myelin_gpas_mode = 'variable'; % sets to 0.02 uF/cm2
    mout.include_rpa = 0; % includes periaxonal resistivity
    mout.myelin_layers = 1; 
end 
if in.print_output
    fprintf('Using params for nrn_model_ver: %s\n',nrn_model_ver)
end
end
% outputs true if cell is PC from L5 or L6, false otherwise
function t_or_f = isL5_6PC(m)
    if isfield(m,'cell_model_name') && ~isempty(m.cell_model_name)
        cell_name = m.cell_model_name;
    else
        cell_name = cellModelNames(m.cell_id,'cell_models_file',m.cell_models_file);
    end
    if m.cell_id > 0 %% BB cell naming convention
        [layer,m_type] = cellNameParser(cell_name);
        if (strcmp(layer,'L5') || strcmp(layer,'L6')) && any(regexp(m_type,'PC','ONCE'))
            t_or_f = true;
        else
            t_or_f = false;
        end
    else % Other
        if ~isempty(regexp(m.cell_model_name,'L5','ONCE')) && ~isempty(regexp(m.cell_model_name,'PC','ONCE'))
           t_or_f = true; 
        else
            t_or_f = false;
        end
    end
    
end

function t_or_f = isPC(m)
    if isfield(m,'cell_model_name') && ~isempty(m.cell_model_name)
        cell_name = m.cell_model_name;
    else
        cell_name = cellModelNames(m.cell_id,'cell_models_file',m.cell_models_file);
    end
    if m.cell_id > 0 % BB cell naming convention
        [~,m_type] = cellNameParser(cell_name);
    else
        m_type = cell_name; % search whole cell name below
    end
    if any(regexp(m_type,'PC','ONCE'))
       t_or_f = true; 
    else
       t_or_f = false;
    end
end
    