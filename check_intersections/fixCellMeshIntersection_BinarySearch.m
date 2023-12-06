function [cell_origin_new, phi_new, fixed] = fixCellMeshIntersection_BinarySearch(cell_origin,cell_normal,...
                                                                                  C,phi,parent_inds,outer_mesh,...                                                                                  
                                                                                  varargin)
%% FIXCELLMESHINTERSECTION_BINARYSEARCH Uses binary search to find minimal translation/rotation
% to move cell within mesh boundaries


in.rel_distance_min = 0.001; % start the minimum as 0.1% of the maximum distance 
in.search_tol = 0.05; % terminate search when relative difference is below this value
in.phi_step = 10; % step for azimuthal rotation search
in.max_iterations = 500; % max number of iterations before terminating search
in.range = [];  % find the range of x,y and z based on the cell self coordinates
in.epsilon = 1e-3; % act as a safety factor
in = sl.in.processVarargin(in,varargin); 
shift_direction = -1 * cell_normal; % assume translation direction is opposite cell_normal 

% set up the boundary for the binary search
if(isempty(in.range) || in.range == 0)
    search_range = range(C(:,3));
else
    search_range = in.range;
end


distance_max = search_range; % get the range of z axis to set up the boundary by height of neuron (max distance in z direction)
distance_min = in.rel_distance_min * distance_max; 

% Check initial position
Cpos = placeCell(cell_origin,cell_normal,C,phi);
[found_intersection,~,~] = checkCellMeshIntersection(outer_mesh,Cpos,parent_inds,...
                                                                 [],0);
if found_intersection    
    fixed = 0;
else    
    fixed = 1; 
    cell_origin_new = cell_origin;
    phi_new = phi; 
    return; 
end
% if there is more than 1% difference for the maximum and minimum
iteration = 0;
phi_new = phi; % initialize with starting phi
while iteration <= in.max_iterations && (distance_max - distance_min)/distance_max > in.search_tol 
    % keep looping until the deviation is below 1 percent
    iteration = iteration + 1;
    current_shift_distance = 0.5 * (distance_max + distance_min); % Assign the current distance as mean of max and min. 
    shift_vector = current_shift_distance * shift_direction; % moving_vector is the moving distance multiplies the moving direction
    cell_origin_new = cell_origin + shift_vector; % update the new_origin by moving the original cell_origin along the moving direction 
    
    % Check azimuthal rotations for this shift distance
    for rotate_phi = 0:in.phi_step:360-in.phi_step
        % conversion so that the angle stays within 360.
        test_phi = mod((rotate_phi + phi), 360);
        % try to place the cell to the new position and detect the
        % intersections
        Cpos = placeCell(cell_origin_new,cell_normal,C,test_phi);
        
        [found_intersection,~,~] = checkCellMeshIntersection(outer_mesh,Cpos,parent_inds,...
            [],0);
        
        % If no intersection is detected, then update the distance_max to
        % lower the upper boundary to the current_distance
        % If intersection is detected, then the distance_min need to
        % increase to the current_distance
        if found_intersection == 0
            fixed = 1;
            phi_new = test_phi;
            %                 fprintf('%g phi worked, fixed = %g, current_shift_distance = %f\n',...
            %                          phi_new,fixed,current_shift_distance);
            break;
        else
            fixed = 0;
        end
    end
    % If for all rotation, there is an intersection, then distance_min need to
    % be updated and set as the current_distance
    if ~fixed
        distance_min = current_shift_distance; % try larger shift distance
    else
        distance_max = current_shift_distance; % try smaller shift distance
    end
    if iteration > in.max_iterations
       fprintf('WARNING: Reached max iterations %g exiting search\n',in.max_iterations); 
    end
end

% phi_new set above

if found_intersection
    fixed = 0;
    fprintf('WARNING: Intersection still found for output position\n');
else
    fixed = 1;
end

end