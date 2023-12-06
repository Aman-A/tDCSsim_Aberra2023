function [cell_data,inds] = extractCellSubMorph(cell_data,comp_type,varargin)
%EXTRACTCELLSUBMORPH Extracts subregion of cell morphology for plotting 
%  
%   Inputs 
%   ------ 
%   cell_data : struct
%               output of loadCellData, includes cell morphology data
%               must include fields: C, comp_types, and parent_inds
%   comp_type : string
%               compartment type, listed in getCellCompInds:
%               'axon','dendrite','apical','basal','soma',or 'all' 
%               also: 'soma_axon', 'soma_dendrite'
%   sectype : string
%               section type, listed in getCellCompInds:
%               'terminal','intermediate','bifurcation','end_points_only', or 'all'
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   cell_data : struct
%               cell_data with only subregion of coordinates in inds and 
%               rest of morphology excluded 
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
in.sectype = 'all';
in = sl.in.processVarargin(in,varargin); 
% cell_data should have comp_types and sectypes
if ~any(strcmp(comp_type,{'all','soma_axon','soma_dendrite'}))
   error(['Must include soma (root section) in morphology at this time',...
          'Use sectype = ''soma_axon'' or ''soma_dendrite''']);  
end
inds = getCellCompInds(cell_data,comp_type,in.sectype);
inds = find(inds); % convert to indices
cell_data.C = cell_data.C(inds,:); 
cell_data.comp_types = cell_data.comp_types(inds); 
cell_data.sectypes = cell_data.sectypes(inds); 
if isfield(cell_data,'secnames')
    cell_data.secnames = cell_data.secnames(inds); 
end
% get connectivity vector (parent_inds)
% parent_inds0 = cell_data.parent_inds(inds(2:end)-1); % assumes 1st index is root
parent_inds0 = cell_data.parent_inds(inds(inds~=1)-1); % soma parent not included in parent_inds
if strcmp(in.sectype,'end_points_only')
%     fprintf('finding end points\n');
%     tic;
    for ii = 1:length(parent_inds0)     
        if (~any(parent_inds0(ii)==inds)) % parent is not a branch/terminal                        
            while true % ascend branch until root is found
                pii = cell_data.parent_inds(parent_inds0(ii)-1); % look at parent's parent
                if (any(pii==inds)) % if reached a comp that is a branch/terminal exit
                    parent_inds0(ii) = pii;
                    break; 
                else
                    parent_inds0(ii) = pii; 
                end
            end        
        end
    end
%     time_elapsed = toc; 
%     fprintf('time elapsed = %f\n',time_elapsed)
end
parent_inds = ones(length(inds)-1,1);
for j = 2:length(parent_inds)
    if (parent_inds0(j) ~= 1) % skip if parent is root (1st index), avoids finding multiple daughter comps
        parent_inds(j) = find(parent_inds0(j) == inds(2:end))+1;
    end
end
cell_data.parent_inds = parent_inds; 