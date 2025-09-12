function pathLoss = computePathloss(positionsUAV, positionsGround, pathLossData)
    % 参数说明：
    % positionsUAV: 无人机位置矩阵，大小为 (3, M, T)
    % positionsGround: 地面发射机位置矩阵，大小为 (2, N, T)
    % pathLossData: 预先计算的路径损耗数据，四维矩阵
    xmin_g = -100;
    ymin_g = -160;
    % 获取尺寸
    [~, M, T] = size(positionsUAV);
    [~, N, ~] = size(positionsGround);
    
    % 初始化路径损耗矩阵
    pathLoss = zeros(M, N, T);
   
    zMin = 10;
    grid_resolution = 5;

    % 遍历所有时隙
    for t = 1:T
        % 获取当前时隙的地面发射机位置
        groundPositions = positionsGround(:, :, t);
        % 获取当前时隙的无人机位置
        uavPositions = positionsUAV(:, :, t);
        
        % 遍历所有无人机和地面发射机
        for m = 1:M
            for n = 1:N
                groundX = groundPositions(1,n);
                groundY = groundPositions(2,n);
                groundXIndex = floor((groundX - xmin_g) / grid_resolution) + 1;
                groundYIndex = floor((groundY - ymin_g) / grid_resolution) + 1;
                groundIndex = (groundYIndex-1)*61+groundXIndex;
            
                % 计算空中位置对应的网格坐标
                airX = uavPositions(1,m);
                airY = uavPositions(2,m);
                airZ = uavPositions(3,m);
                airXIndex = floor((airX - xmin_g) / grid_resolution) + 1;
                airYIndex = floor((airY - ymin_g) / grid_resolution) + 1;
                airZIndex = floor((airZ - zMin) / grid_resolution)+1;
            
                % 检查索引是否在有效范围内
                if groundXIndex < 1 || groundXIndex > size(pathLossData, 1) || ...
                   groundYIndex < 1 || groundYIndex > size(pathLossData, 2) || ...
                   airXIndex < 1 || airXIndex > size(pathLossData, 1) || ...
                   airYIndex < 1 || airYIndex > size(pathLossData, 2) || ...
                   airZIndex < 1 || airZIndex > size(pathLossData, 3)
                   pathLoss(m,n,t) = 500;
                   return
                end
            
                % 获取地面位置的路径损耗值
                pathLoss(m,n,t) = pathLossData(airXIndex, airYIndex, airZIndex, groundIndex);
            end
        end
    end
end

