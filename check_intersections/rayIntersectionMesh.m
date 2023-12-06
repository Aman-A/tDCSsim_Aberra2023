function [intersect,t] = rayIntersectionMesh(origin, direction, mesh)
%% This function checks whether the ray is intersecting with a mesh
% Input: origin, direction, mesh
% Output: flag, t
% flag is 1 when the ray is intersecting with the mesh
% t is the nonegative distance from the intersection point to the direction

layer_faces = mesh.faces;
layer_vertices = mesh.vertices;
distance = [];
intersect = 0;
for face_index = 1:size(layer_faces,1)
    temp_face = layer_faces(face_index,:);
    temp_vertices = layer_vertices(temp_face,:);
    p0 = temp_vertices(1,:);
    p1 = temp_vertices(2,:);
    p2 = temp_vertices(3,:);
    [~, u, v, t] = fastRayTriangleIntersection(origin, direction, p0, p1, p2);
    if u > 0 && v > 0 && u + v < 1
        distance = [distance, t];
        if (t > 0)
            intersect = 1;
        end
    end
end
t = min(distance(distance > 0));

end