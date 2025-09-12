function positions = generate_path_points(start_position, end_position, num_points)
    % start_position 和 end_position 是 3x1 向量，表示 UAV 的起点和终点
    % num_points 是需要生成的路径点数量

    % 初始化位置数组
    positions = zeros(2, num_points);

    % 计算每个维度的增量
    delta = (end_position - start_position) / (num_points - 1);

    % 生成路径点
    for i = 1:num_points
        positions(:, i) = start_position + delta * (i - 1);
    end
end
