function cell_depth = getCellDepth(cell_origin,gm,wm)
%% This is a function for finding the cell depth
% This is achieved by finding the 1 nearest neighbor of the cell origin to 
% the white matter and gray matter
cell_depths = zeros(size(cell_origin,1),1);

for i = 1:size(cell_origin,1)
    Idx = knnsearch(cell_origin, gm.vertices);

    gm_nearest_idx = find(Idx == 1);
    gm_nearest = gm.vertices(gm_nearest_idx);
    gm_distance = vmag(gm_nearest - cell_origin);

    Idx = knnsearch(cell_origin, wm.vertices);

    wm_nearest_idx = find(Idx == 1);
    wm_nearest = wm.vertices(wm_nearest_idx);
    wm_distance = vmag(wm_nearest - cell_origin);
    cell_depth = gm_distance/(gm_distance + wm_distance);
end

end