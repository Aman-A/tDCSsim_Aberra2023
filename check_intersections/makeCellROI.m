function cellROI = makeCellROI(C,ROI_expand_factor)
% Makes box around extremes of cell morphology coordinates, expanded by
% ROI_expand_factor
minimum = min(C,[],1);
maximum = max(C,[],1);

% [xmin1 xmax1 ymin1 ymax1 zmin1 zmax1]
cellROI = [minimum(1) maximum(1) minimum(2) maximum(2) minimum(3) maximum(3)]; 
% Scale by expand_factor
% cellROI = scale(cellROI,expand);
for i = [1,3,5] % expand ROI to include points just outside
    cellROI(i:i+1) = cellROI(i:i+1)+ROI_expand_factor*diff(cellROI(i:i+1))*[-1,1];
end
end