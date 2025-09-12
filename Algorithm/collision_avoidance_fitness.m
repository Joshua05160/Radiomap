function fitness = collision_avoidance_fitness(path_points, building_heights_and_floors)
    % path_points: 3 x N x T matrix, where each column represents a point [x; y; z]
    % building_heights_and_floors: Mx3 matrix, where each row represents [building_id, height, floors]

    % Initialize fitness
    fitness = 0;
    load("building_data.mat");
    building_ids = unique(filtered_data.BuildingID);
    num_buildings = length(building_ids);

    for t = 1:size(path_points, 3)
        % Loop through each path point
        for i = 1:size(path_points, 2)
            point = path_points(:, i, t);
            x = point(1);
            y = point(2);
            z = point(3);
            if t < size(path_points, 3)
                next_point = path_points(:, i, t+1);
                mid_point = (point + next_point) / 2;
                x_mid = mid_point(1);
                y_mid = mid_point(2);
                z_mid = mid_point(3);
            end
            % Check each building
            for b = 1:length(building_ids)
                building_id = building_ids(b);
                % Get the building vertices
                building_vertices = filtered_data(filtered_data.BuildingID == building_ids(b), 2 :3);
                building_vertices = building_vertices(1:end-1,:);
                % Get the building height
                building_info = building_heights_and_floors(building_heights_and_floors(:, 1) == building_ids(b), :);
                building_height = building_info(2);
                % Check if the point is inside the building footprint
                [in1, on] = inpolygon(x, y, building_vertices.X, building_vertices.Y);
                in2 = 0;
                if t < size(path_points, 3)
                    [in2, on] = inpolygon(x_mid, y_mid, building_vertices.X, building_vertices.Y);
                end
                
                if in1
                    % Check for collision
                    height_diff = z - building_height;
                    if height_diff < 0
                        fitness = fitness + height_diff;
                    end
                end
                if in2
                   height_diff = z_mid - building_height;
                    if height_diff < 0
                        fitness = fitness + height_diff;
                    end
                end
        end
    end
    end
    fitness = fitness*10;
end
