function [varargout] = removeMeshFaces(vertices,faces,inds_remove,varargin)
%REMOVEMESHFACES Removes mesh faces for input indices and outputs mesh
%  
%   Inputs 
%   ------ 
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
in.vert_inds = []; % vertices to remove after reordering using faces
in = sl.in.processVarargin(in,varargin); 
f = faces(~inds_remove,:);
if size(vertices,2) == 1 % convert Nx1 vector to Nx3 to match vertices arrays
    p1 = repmat(vertices(f(:,1),:),1,3);
    p2 = repmat(vertices(f(:,2),:),1,3);
    p3 = repmat(vertices(f(:,3),:),1,3);
else
    p1 = vertices(f(:,1),:); 
    p2 = vertices(f(:,2),:); 
    p3 = vertices(f(:,3),:);
end
p_all = [p1';p2';p3'];
p_all = reshape(p_all,3,length(p1)*3)'; % [f1 pt1;f1 pt2;f1 pt3; f2 pt1; f2 pt2;...]
if isempty(in.vert_inds) % remove vertices for corresponding faces
    [p_un,vert_inds,ic] = unique(p_all,'rows','stable');
    facesN = reshape((1:size(p_all,1)),3,length(p1))'; % new triangle faces
    facesN = ic(facesN); % renumber for unique vertices
else % remove vertices based on input vertex indices
    vert_inds = in.vert_inds; 
    p_un = p_all(vert_inds,:); 
    if all(~diff(p_un,[],2),'all')
       p_un = p_un(:,1); % take first column if all identical 
    end
end
verticesN = p_un; % retain only unique points
if nargout == 3
    varargout = {verticesN,facesN,vert_inds};
elseif nargout == 2
    varargout = {verticesN,facesN};
elseif nargout == 1
    varargout = {verticesN}; 
end
end