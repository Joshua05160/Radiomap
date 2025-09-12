function [positions, velocities] = initialize_particles_multi(start_position, end_position, M, T, numParticles, lastBestPosition, speedLimit)
    % start_position 和 end_position 是 3xM 矩阵，表示 M 架 UAV 的起点和终点
    % M 是无人机的数量
    % T 是时间步长
    % numParticles 是粒子的数量
    
    % 初始化粒子位置
    positions = zeros(3, M, T, numParticles);
    
    % 设置第一个粒子的位置为上次优化的最优位置
    positions(:, :, :, 1) = lastBestPosition;
    
    for m = 1:M
        % 设置所有粒子的初始位置为相同的起点和终点
        for i = 1:numParticles
            positions(:, m, 1, i) = start_position(:, m); % 初始化起点
            positions(:, m, T, i) = end_position(:, m);   % 初始化终点
        end
        
        % 计算起点和终点之间的连线向量
        line_vector = end_position(:, m) - start_position(:, m);
        line_length = norm(line_vector);
        unit_vector = line_vector / line_length;
        
        % 在连线上等距离选择三个垂直于该连线的面
        for i = 2:numParticles
            % 初始化随机点的中间位置
            random_points = zeros(3, 3); % 存储三个随机点
            for j = 1:3
                % 计算垂直平面的中心点
                plane_center = start_position(:, m) + j * line_length / 4 * unit_vector;
                
                % 随机选择垂直于连线的面的两个正交向量
                if unit_vector(1) ~= 0 || unit_vector(2) ~= 0
                    orthogonal_vector1 = [-unit_vector(2); unit_vector(1); 0];
                else
                    orthogonal_vector1 = [1; 0; 0];
                end
                orthogonal_vector2 = cross(unit_vector, orthogonal_vector1);
                
                % 在垂直平面上随机选择一个点
                random_point = plane_center + (rand - 0.5) * 4 * speedLimit * orthogonal_vector1 + (rand - 0.5) * 4 * speedLimit * orthogonal_vector2;
                random_points(:, j) = random_point;
            end
            
            % 通过插值将起点、终点和三个随机点连接起来并等距离采样
            sample_points = [start_position(:, m), random_points, end_position(:, m)];
            for t = 1:T
                interpolated_position = interp1(linspace(0, 1, 5), sample_points', (t - 1) / (T - 1))';
                positions(:, m, t, i) = interpolated_position;
            end
        end
    end
    
    % 初始化粒子速度
    velocities = 2 * 10 * rand(3, M, T, numParticles) - 10;  % 在[-10, 10]范围内生成随机速度
    
    % 固定起点和终点的速度为0
    velocities(:, :, 1, :) = 0;
    velocities(:, :, T, :) = 0;
end
