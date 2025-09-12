function smoothed_positions = smooth_trajectory(positionsUAV, window_size)
    % 输入:
    % positionsUAV: 原始轨迹数据，尺寸为 [维度, M, T]，维度可以是2或3
    % window_size: 移动平均窗口大小，必须是奇数

    % 输出:
    % smoothed_positions: 平滑后的轨迹数据，尺寸与 positionsUAV 相同

    % 获取轨迹数据的尺寸
    [dim, M, T] = size(positionsUAV);
    
    % 预分配空间
    smoothed_positions = zeros(dim, M, T);
    
    % 对每条轨迹进行平滑处理
    for m = 1:M
        for d = 1:dim
            % 获取当前维度的轨迹数据
            trajectory = squeeze(positionsUAV(d, m, :));
            smoothed_trajectory = trajectory;
            % 使用移动平均法进行平滑处理
            smoothed_trajectory(2:T-1) = movmean(trajectory(2:T-1), window_size);
            
            % 存储平滑后的轨迹数据
            smoothed_positions(d, m, :) = smoothed_trajectory;
        end
    end
end
