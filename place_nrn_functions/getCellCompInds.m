function inds = getCellCompInds(cell_data,comp_type,sectype)
if nargin == 0
   cell_data = loadCellData(6,'maxH',{'comp_types','sectypes'}); 
   comp_type = 'axon';
   sectype = 'all';
end
sectypes = cell_data.sectypes; 
comp_types = cell_data.comp_types;
% filter by compartment type:
% 0 - soma, 1 - axon, 2 - node, 3 - myelin, 4 - unmyelin, 5 - basal, 6 - apic
if strcmp(comp_type,'axon')
    comp_type_inds = comp_types > 0 & comp_types < 5;
elseif strcmp(comp_type,'dendrite')
    comp_type_inds = comp_types == 5 | comp_types == 6;
elseif strcmp(comp_type,'apical')
    comp_type_inds = comp_types == 6;
elseif strcmp(comp_type,'basal')
    comp_type_inds = comp_types == 5;
elseif strcmp(comp_type,'soma')
    comp_type_inds = comp_types == 0;
elseif strcmp(comp_type,'soma_axon')
    comp_type_inds = comp_types >= 0 & comp_types < 5;
elseif strcmp(comp_type,'soma_dendrite')
    comp_type_inds = comp_types == 0 | comp_types == 5 | comp_types == 6;
elseif strcmp(comp_type,'all') || strcmp(comp_type,'any')
    comp_type_inds = true(size(comp_types));
else
    error('comp_type ''%s'' not defined\n',comp_type); 
end
% filter by section type:
% 1 = soma (no parent)
% 2 = termination from non-bifurcation (parent is a 1,3, or 6) OR far from bifurcation (parent is 4 or 7, L > length constant) ** contains termination point (1) ** 
% 3 = intermediate from non-bifurcation (parent is 1,3, or 6)
% 4 = parent side of bifurcation (child is a 5,6, or 7) (parent is 3 OR far from bifurcation (4 or 7)) ** contains parent side bifurcation point (1)** 
% 5 = child termination from bifurcation (parent is a 1 or 4 and L < length constant) ** contains termination point (1)** 
% 6 = child side of bifurcation (parent is a 1 or 4) **contains child side bifurcation point (0)**
% 7 = parent bifurcation and child of bifurcation (parent is a 1 or 4, child is a 5,6, or 7) ** contains parent side bifurcation point (1) **
if strcmp(sectype,'terminal')
   sectype_inds =  sectypes == 2 | sectypes == 5;
elseif strcmp(sectype,'intermediate')
   sectype_inds =  sectypes == 3;
elseif strcmp(sectype,'bifurcation')
    sectype_inds =  sectypes == 4 | sectypes > 5;
elseif strcmp(sectype,'end_points_only')
    sectype_inds = sectypes == 1 | sectypes == 2 | sectypes == 5 | ... % terminals
                    sectypes == 4 | sectypes == 7; % parent bifurcation or parent and child bifurcation    
elseif strcmp(sectype,'all') || strcmp(sectype,'any')
    sectype_inds = true(size(sectypes));
else
    error('sectype ''%s'' not defined\n',sectype); 
end
inds = comp_type_inds & sectype_inds;

end