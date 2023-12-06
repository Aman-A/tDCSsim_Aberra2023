function p = plotDefPatch(patch_struct)
% PLOTDEFPATCH input patch structure, plots patch with default settings to current axis
% patch_struct:
%   faces: [N x 3 double]
%   vertices: [N x 3 double]
% TODO allow patch arguments as function args using varargin
def_edge_color = 'none';
def_face_color = [0.8 0.8 0.8];
def_alpha = 1; 
def_lighting = 'gouraud';
% copy into new struct for plotting (avoids conflict when patch_struct has
% extraneous fields)
if isfield(patch_struct,'faces')
   patch_struct_plot.faces = patch_struct.faces;
elseif isfield(patch_struct,'triangles')
    patch_struct_plot.faces = patch_struct.triangles;
elseif isfield(patch_struct,'Faces')
    patch_struct_plot.faces = patch_struct.Faces;
end
if isfield(patch_struct,'vertices')
    patch_struct_plot.vertices = patch_struct.vertices;
elseif isfield(patch_struct,'points')
    patch_struct_plot.vertices = patch_struct.points;
elseif isfield(patch_struct,'nodes')
    patch_struct_plot.vertices = patch_struct.nodes;
elseif isfield(patch_struct,'Vertices')
    patch_struct_plot.vertices = patch_struct.Vertices;
end
if isa(patch_struct,'triangulation')
    patch_struct_plot.faces = patch_struct.ConnectivityList;
    patch_struct_plot.vertices = patch_struct.Points;
end
p = patch(patch_struct_plot,'EdgeColor',def_edge_color,'FaceColor',def_face_color,...
    'FaceAlpha',def_alpha,'FaceLighting',def_lighting);
axis equal; 
ax = gca;
types = arrayfun(@(x) x.Type,ax.Children,'UniformOutput',0);
if ~any(strcmp(types,'light'))
    camlight head; 
end