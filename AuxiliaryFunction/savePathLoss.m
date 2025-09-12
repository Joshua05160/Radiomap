% 定义文件夹路径

% 定义网格分辨率
% 定义网格分辨率
grid_resolution = 5;

% 生成x和y范围
x_range = -100:grid_resolution:200;
y_range = -160:grid_resolution:160;

% 生成网格坐标并交换参数位置
[Y, X] = meshgrid(y_range, x_range);

% 将网格坐标展开为向量并组合
coordinates = [X(:), Y(:)];
    % 定义 x, y, z 的范围和分辨率
    xMin = -97.5;
    xMax = 197.5;
    yMin = -157.5;
    yMax = 157.5;
    zMin = 10;
    zMax = 55;
    resolution = 5;
    
    % 计算 x, y, z 方向的大小
    X = round((xMax - xMin) / resolution) + 1;
    Y = round((yMax - yMin) / resolution) + 1;
    Z = round((zMax - zMin) / resolution) + 1;
    M = 3965; % 代表一个地面发射机的坐标
    
    % 初始化四维数组
    pathLossData = zeros(X, Y, Z, M);

for m=1:M
    folderPath = sprintf('（%d，%d）',coordinates(m,1),coordinates(m,2)); % 文件夹路径
    % 获取文件夹中所有文件的列表
    filePattern = fullfile(folderPath, '*.txt'); % 假设文件是以.txt结尾
    files = dir(filePattern);
    % 读取文件并存储数据
    for k = 1:length(files)
        baseFileName = files(k).name;
        fullFileName = fullfile(files(k).folder, baseFileName);
        
        % 读取文件内容，跳过前11行头信息
        data = readmatrix(fullFileName);
       
        % 提取路径损耗值并存储到四维数组中
        for i = 1:size(data, 1)
            x = round((data(i, 1) - xMin) / resolution) + 1;
            y = round((data(i, 2) - yMin) / resolution) + 1;
            z = round((data(i, 3) - zMin) / resolution) + 1;
            % 确保坐标在有效范围内
            pathLoss = data(i, 6);
            pathLossData(x, y, z, m) = pathLoss;
        end
    end
end
% 保存数据到 MAT 文件
save('pathLossData.mat', 'pathLossData');
