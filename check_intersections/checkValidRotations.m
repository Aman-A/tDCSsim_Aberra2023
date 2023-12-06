function [valid_rotation,rotate_degrees] = checkValidRotations(cell_origin,...
                                            cell_normal,C,phi,parent_inds,...
                                            outer_mesh,varargin)

in.phi_step = 30;      
in = sl.in.processVarargin(in,varargin); 

valid_rotation = zeros(1,360/in.phi_step); % the boolean array specify if the 12 rotations are ok
rotate_degrees = zeros(1,360/in.phi_step); % the actual angle that is checked for the 12 rotations
rotate_count = 1;

for rotate_phi = 0:in.phi_step:360-in.phi_step
        % conversion so that the angle stays within 360.
        test_phi = mod((rotate_phi + phi), 360);
        rotate_degrees(rotate_count) = test_phi;
        % try to place the cell to the new position and detect the
        % intersections
        Cpos = placeCell(cell_origin,cell_normal,C,test_phi);
        
        [found_intersection,~,~] = checkCellMeshIntersection(outer_mesh,Cpos,...
                                                            parent_inds,...
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
            valid_rotation(rotate_count) = 1;
        else
            fixed = 0;
            valid_rotation(rotate_count) = 0;
            
        end
        rotate_count = rotate_count + 1;
end

end