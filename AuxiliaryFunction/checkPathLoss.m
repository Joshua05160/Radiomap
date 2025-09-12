% 当前文件夹路径
currentFolderPath = 'PATH_loss';

% 获取当前文件夹下的所有文件夹列表
currentFolders = dir(currentFolderPath);
currentFolders = currentFolders([currentFolders.isdir]); % 只保留文件夹
currentFolderNames = {currentFolders.name};
currentFolderNames = currentFolderNames(~ismember(currentFolderNames, {'.', '..'}))% 排除 '.' 和 '..'

% 给定的文件夹列表（假设给定列表中的文件夹名称格式为 '(x,y)'）
givenFolders = cell(3965, 1);
for i = 1:3965
    givenFolders{i} = sprintf('（%d，%d）', coordinates(i, 1), coordinates(i, 2));
end
% 初始化缺失文件夹矩阵
missingFolders = {};

% 比较文件夹，找出缺失的文件夹
for i = 1:length(givenFolders)
    if ~ismember(givenFolders{i}, currentFolderNames)
        missingFolders{end+1} = givenFolders{i};
    end
end

% 将缺失文件夹转换为矩阵形式
missingFoldersMatrix = char(missingFolders);

% 显示缺失文件夹矩阵
disp('缺失的文件夹:');

