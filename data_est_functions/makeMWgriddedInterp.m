% Input either single vector of theta,phi, and thresholdE values for cell
% or cell array of thetas,phis, and thresholdE vectors for multiple cells
% outputs griddedInterpolant object or objects in cell array, respectively
% any threshold value can then be extracted by sampling the interpolated
% grid at a given (theta,phi):
% THRESH_ints{cell_id}(theta,phi)
% or THRESH_ints(theta,phi)
function DATA_ints = makeMWgriddedInterp(thetas,phis,data,varargin)
in.interp_mode = 'linear';
in = sl.in.processVarargin(in,varargin);
if iscell(data)
    num_cells = length(data); % input single or multiple cells' threshold vectors in cell array
    if ~iscell(thetas) % thetas/phis are vectors, not cell arrays
        thetas = repmat({thetas},num_cells,1);
        phis = repmat({phis},num_cells,1);
    end
else
    num_cells = 1; % Input a single cell's threshold vector
    thetas = {thetas}; phis = {phis}; data = {data}; % for accessing in for loop below
end
DATA_ints = cell(num_cells,1); % initialize cell array for threshold interp grids
for i = 1:num_cells
    if(isempty(thetas{i}) && isempty(phis{i}) && isempty(data{i}))
        DATA_ints{i} = [];
    else
        theta_i = thetas{i}; phi_i = phis{i}; threshEs_i = data{i};
        num_lat = length(unique(theta_i)); num_long = length(unique(phi_i))+1;
        
        % create grids for interpolation
        
        THETA=zeros(num_lat,num_long);
        PHI=zeros(num_lat,num_long);
        DATA=zeros(num_lat,num_long);
        
        THETA(2:num_lat-1,1:num_long-1)=reshape(theta_i(2:end-1),num_long-1,num_lat-2)';
        THETA(2:num_lat-1,num_long)=THETA(2:num_lat-1,1);
        THETA(1,:)=theta_i(1);
        THETA(num_lat,:)=theta_i(end);
        
        PHI(2:num_lat-1,1:num_long-1)=reshape(phi_i(2:end-1),num_long-1,num_lat-2)';
        PHI(1,:)=PHI(2,:);
        PHI(num_lat,:)=PHI(num_lat-1,:);
        ind=(PHI(:,1:num_long-1)>=180);
        PHI(ind)=PHI(ind)-360;
        [~,sort_ind]=sort(PHI(1,1:num_long-1));
        PHI=[PHI(:,sort_ind),PHI(:,sort_ind(1))+360-eps];
        
        DATA(2:num_lat-1,1:num_long-1)=reshape(threshEs_i(2:end-1),num_long-1,num_lat-2)';
        DATA(1,:)=threshEs_i(1);
        DATA(num_lat,:)=threshEs_i(end);
        DATA=[DATA(:,sort_ind),DATA(:,sort_ind(1))];
        
        THETA=[THETA(:,sort_ind),THETA(:,sort_ind(1))];
        
        THRESH_int_i = griddedInterpolant(THETA,PHI,DATA,in.interp_mode); % create gridded interpolant object for ith cell
        
        DATA_ints{i} = THRESH_int_i;
    end
end

end