function [found_intersection,mesh_faces,mesh_vertices] = checkCellMeshIntersection(outer_mesh,Cpos,parent_inds,...
                                                           cellROI_expand_factor,crop_outer_mesh)
% CHECKCELLMESHINTESRECTION Checks if cell at single position intersects
% with nearby portion of surface mesh (e.g. GrayMatter)
if nargin < 5
    crop_outer_mesh = 1; 
end
if crop_outer_mesh
    [mesh_faces,mesh_vertices] = getMeshNearCell(outer_mesh,Cpos,cellROI_expand_factor); 
else % skip cropping (already cropped externally)
    mesh_faces = outer_mesh.faces; 
    mesh_vertices = outer_mesh.vertices; 
end
found_intersection = 0;
for i = 1:size(mesh_faces,1) % loop through cropped GM mesh faces
    % vertices of face i
    v1 = mesh_vertices(mesh_faces(i,1),:);
    v2 = mesh_vertices(mesh_faces(i,2),:);
    v3 = mesh_vertices(mesh_faces(i,3),:);
    for j = 2:size(Cpos,1) % loop through line segments of cell morphology
        p1 = Cpos(parent_inds(j-1),:); % parent of coord j
        p2 = Cpos(j,:); % coord j
        [intersect,~,~,~] = fastRayTriangleIntersection(p1,p2-p1,v1,v2,v3); 
        if(intersect == 1)
            found_intersection = 1;
            return; % break out of loop for this mesh face
        end
    end
end
end