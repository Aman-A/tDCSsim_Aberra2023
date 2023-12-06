function [lower_intersect, upper_intersect,lower_distance, upper_distance] = ...
        checkInsideOriginalLayer(Cposi,cell_origin,cell_normal,upper_boundary_surf,...
                                lower_boundary_surf,varargin)
    %% Input
    % -- cell_origin: 1 by 3 vector with coordinates of (x,y,z)
    % -- cell_normal: 1 by 3 vector with vector component of (x,y,z)
    % -- partial_lb: the mesh of the lower boundary mesh
    % Output
    % -- lower_intersect, upper_intersect: whether the cell normal intersects with the lower boundary
    % or the inverted cell normal intersects with the 
    % -- distance_to_layer: the distance to the lower boundary 
    % Shoot an inverted normal and see if it intersects with the layer
    % boundary mesh
     % Check if the cell normal intersect the lower boundary surface
    in.cellROI_expand_factor = 0.5;
    in.min_nearest_faces = 10;
    in = sl.in.processVarargin(in,varargin);
    
    cellROI_expand_factor = in.cellROI_expand_factor;
    min_nearest_faces = in.min_nearest_faces;
    
    [lower_layer_faces,lower_layer_vertices] = getMeshNearCell(lower_boundary_surf,...
                                                Cposi,cellROI_expand_factor,...
                                                'min_nearby_faces',min_nearest_faces); 
    partial_lb = struct('faces',lower_layer_faces,'vertices',lower_layer_vertices);  
    [lower_intersect,lower_distance] = rayIntersectionMesh(cell_origin, cell_normal, partial_lb);
    
    % Check if the inverted cell normal intersect the upper boundary surface 
    [upper_layer_faces,upper_layer_vertices] = getMeshNearCell(upper_boundary_surf,...
                                                Cposi,cellROI_expand_factor,...
                                                'min_nearby_faces',min_nearest_faces); 
    partial_ub = struct('faces',upper_layer_faces,'vertices',upper_layer_vertices);  
    [upper_intersect,upper_distance] = rayIntersectionMesh(cell_origin, -1*cell_normal, partial_ub);
    
end