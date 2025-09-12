function new_ALPHA = interpolate_trajectory(ALPHA, new_T)
    % ALPHA: 输入的轨迹矩阵，包含原始采样点
    % new_T: 新的采样点数
    
    % 初始化新的 ALPHA 矩阵
    Msize = size(ALPHA,1);
    % 对每一段轨迹进行插值
    if Msize == 3
        new_ALPHA = zeros(3, size(ALPHA, 2), new_T);
        for i = 1:size(ALPHA, 2)
            % 获取当前轨迹段的 x 和 y 坐标
            x = ALPHA(1, i, :);
            y = ALPHA(2, i, :);
            z = ALPHA(3, i, :);
 
            % 将坐标转换为列向量
            x = squeeze(x);
            y = squeeze(y);
            z= squeeze(z);
            % 使用线性插值将采样点数从 T 增加到 new_T
            new_x = interp1(1:numel(x), x, linspace(1, numel(x), new_T));
            new_y = interp1(1:numel(y), y, linspace(1, numel(y), new_T));
            new_z = interp1(1:numel(z), z, linspace(1, numel(z), new_T));

            % 更新新的 ALPHA 矩阵
            new_ALPHA(:, i, :) = [new_x; new_y;new_z];
        end
    end
    
    
    if Msize ==2
            new_ALPHA = zeros(2, size(ALPHA, 2), new_T);

        for i = 1:size(ALPHA, 2)
            % 获取当前轨迹段的 x 和 y 坐标
            x = ALPHA(1, i, :);
            y = ALPHA(2, i, :);
           % z = ALPHA(3, i, :);
 
            % 将坐标转换为列向量
            x = squeeze(x);
            y = squeeze(y);
            %z= squeeze(z);
            % 使用线性插值将采样点数从 T 增加到 new_T
            new_x = interp1(1:numel(x), x, linspace(1, numel(x), new_T));
            new_y = interp1(1:numel(y), y, linspace(1, numel(y), new_T));
           % new_z = interp1(1:numel(z), z, linspace(1, numel(z), new_T));

            % 更新新的 ALPHA 矩阵
            new_ALPHA(:, i, :) = [new_x; new_y];
        end
    end
end
