function data_out = processNrnCompData(data,cell_ids,nrn_model_ver,...
                                       func,comp_type,sectype,varargin)
%PROCESSNRNCOMPDATA Process neuron recording data from multiple 
% compartments to result in a single value per simulation (row)
%  
%   Inputs 
%   ------ 
%   data : 1 x num_cells cell array
%          Each cell contains num_simulations x num_compartments matrix
%          of recorded data values, OR each cell can be a 1 x num_pops cell
%          array in which each element contains the num_simulations x
%          num_compartments matrices
%   cell_ids : 1 x num_cells vector
%              cell ids (cell_model_names = getCellModelNames(cell_ids))
%   nrn_model_ver : string
%                   neuron model version (see outputMorphParams)
%   func : string or numeric
%          Function to apply to compartment data from each simulation,
%          either: 'max', 'maxabs', 'mean', 'median','min','maxabs_ind'.
%          'maxabs_ind' extracts value from each simulation at the single 
%          compartment with the maxabs value across all simulations. If
%          numeric and < 1, calculates quantile of func. Input 'none' if
%          you just want to extract data in comp_type/sectype compartments
%          without applying a function to collapse data across compartments
%   comp_type : string
%               compartment type to include in func computation, see
%               getCellCompInds for definitions
%   sectype : string
%             section type to include in func computation, see
%             getCellCompInds for definitions
%   Optional Inputs 
%   --------------- 
%   cell_models_file : string (default = 'cell_models')
%                      *.mat file in mat/cell_data to load table of
%                      cell_model_names, used in cellModelNames
%   Outputs 
%   ------- 
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
if nargin < 3
   func = 'max';
   comp_type = 'axon';
   sectype = 'terminal';
end
in.cell_models_file = 'cell_models';
in.print_level = 0; 
in.cell_data = {}; % cell_data corresponding to cell_ids vector
in = sl.in.processVarargin(in,varargin); 
num_cells = length(cell_ids);
data_out = cell(size(data)); 
if in.print_level > 0
    fprintf('Processing data from %g cells with func: ''%s'', comp_type: ''%s'', sectype: ''%s''\n',...
            num_cells,func,comp_type,sectype); 
end
for i = 1:num_cells
    if iscell(data{i})
        num_recs = size(data{i}{1},2); % multiple rotations per cell, assume all same number of recording locations  
    else
        num_recs = size(data{i},2); % single rotation per cell
    end
    if num_recs > 1 % only run if more than one value recorded per simulation
        if isempty(in.cell_data)
            cell_datai = loadCellData(cell_ids(i),nrn_model_ver,...
                                        {'comp_types','sectypes','secnames'},...
                                        'cell_models_file',in.cell_models_file); 
        else
            cell_datai = in.cell_data{i}; 
        end
        % get desired compartment indices (binary vector)
        inds = getCellCompInds(cell_datai,comp_type,sectype);         
        if any(inds)
            if iscell(data{i}) % for multiple nrn_pops
                data_out{i} = cellfun(@(x) x(:,inds),data{i},'UniformOutput',0); 
            else
                data_out{i} = data{i}(:,inds); 
            end
        else % skip if no compartments found
            data_out{i} = []; 
        end
        % apply function to data from relevant compartments
        if ~strcmp(func,'none') && ~isempty(data_out{i})
            if iscell(data_out{i})
                data_out{i} = cellfun(@(x) apply_func(x,func),data_out{i},'UniformOutput',0);
            else
                data_out{i} = apply_func(data_out{i},func);
                if strcmp(func,'ind_maxabs')
                    inds_lin = find(inds); 
                    data_out{i} = inds_lin(data_out{i});
                end
            end
        end
    else
        data_out{i} = data{i}; 
    end
end
% if length(data_out) == 1
%    data_out = data_out{1};
% end
end