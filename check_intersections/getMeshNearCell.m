function [mesh_faces, mesh_vertices, cellROI] = getMeshNearCell(outer_mesh,Cpos,cellROI_expand_factor,...
                                                varargin)
in.min_nearby_faces = 20; % min number of faces in cropped outer_mesh
in.cellROI_expand_factor_step = 0.5; % step size to increase if mesh has too few faces
in = sl.in.processVarargin(in,varargin);     
mesh_faces = []; 
while size(mesh_faces,1) < in.min_nearby_faces
    % get box at min/max coordinate in x, y, and z coordinate (global)    
    cellROI = makeCellROI(Cpos,cellROI_expand_factor); % (cropping_vertices)
    % Crop outer mesh with cellROI
    [mesh_vertices,mesh_faces,~,~] = clipMeshVertices(outer_mesh,cellROI);

    cellROI_expand_factor = cellROI_expand_factor + in.cellROI_expand_factor_step;
end
end