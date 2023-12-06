function varargout = boxDelete(Points,Box)
%Delete points (X,Y,Z,additionalparams...) within rectangular region
%defined by Box = [xmin xmax ymin ymax zmin zmax]

xmin = Box(1);
xmax = Box(2);
ymin = Box(3);
ymax = Box(4);
zmin = Box(5);
zmax = Box(6);

% compute indices of points inside Noise box
xOk = Points(:,1) >= xmin & Points(:,1) <= xmax;
yOk = Points(:,2) >= ymin & Points(:,2) <= ymax;
zOk = Points(:,3) >= zmin & Points(:,3) <= zmax;

% keep only points inside box
inds = find(xOk & yOk & zOk);  % indices of points inside Box
ind = setdiff(1:size(Points,1),inds); % indices of points outside Box
Points(xOk & yOk & zOk,:) = []; % remove points inside noise box
% Output variables
varargout{1} = Points;
if nargout == 2
    varargout{2} = ind; 
end
end