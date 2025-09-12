function extractAndVisualizeBuildingsFromOSM(osmFilePath, outputMatFile)
% 从OSM文件中提取建筑物坐标并可视化（转换为米制坐标系统）
% 
% 功能描述:
%   1. 解析OSM XML文件，提取建筑物数据
%   2. 将经纬度坐标转换为相对于左下角原点的米制坐标
%   3. 保存数据为.mat文件
%   4. 生成建筑物俯视图和分布密度热图
%
% 输入参数:
%   osmFilePath: OSM文件路径
%   outputMatFile: 输出的.mat文件路径 (可选，默认为'buildings.mat')
%
% 输出数据结构:
%   buildings(i).id: 建筑物编号
%   buildings(i).type: 建筑物类型
%   buildings(i).latlon: 原始经纬度坐标 [纬度, 经度]
%   buildings(i).x_meters: 东向距离 (米，相对于左下角原点)
%   buildings(i).y_meters: 北向距离 (米，相对于左下角原点)
%   buildings(i).coordinates_meters: [x_meters, y_meters]
%
% 坐标系统说明:
%   - 原点位于所有建筑物的左下角（最小经度，最小纬度）
%   - X轴指向东方，单位为米
%   - Y轴指向北方，单位为米
%   - 使用球面近似转换：1度纬度≈111320米，1度经度≈111320*cos(纬度)米
%
% 使用示例:
%   extractAndVisualizeBuildingsFromOSM('map.osm', 'buildings.mat');
%   extractAndVisualizeBuildingsFromOSM('map.osm');  % 默认输出文件名

if nargin < 2
    outputMatFile = 'buildings.mat';
end

fprintf('正在读取OSM文件: %s\n', osmFilePath);

% 读取OSM XML文件
try
    osmData = xmlread(osmFilePath);
catch
    error('无法读取OSM文件，请检查文件路径');
end

% 获取所有节点
nodes = osmData.getElementsByTagName('node');
nodeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

fprintf('正在解析 %d 个节点...\n', nodes.getLength());

% 创建节点ID到坐标的映射
for i = 0:nodes.getLength()-1
    node = nodes.item(i);
    nodeId = char(node.getAttribute('id'));
    lat = str2double(node.getAttribute('lat'));
    lon = str2double(node.getAttribute('lon'));
    nodeMap(nodeId) = [lat, lon];
end

% 获取所有way元素
ways = osmData.getElementsByTagName('way');
buildings = [];
buildingCount = 0;

fprintf('正在解析 %d 个路径，查找建筑物...\n', ways.getLength());

for i = 0:ways.getLength()-1
    way = ways.item(i);
    
    % 检查是否为建筑物
    tags = way.getElementsByTagName('tag');
    isBuilding = false;
    buildingType = 'unknown';
    
    for j = 0:tags.getLength()-1
        tag = tags.item(j);
        key = char(tag.getAttribute('k'));
        value = char(tag.getAttribute('v'));
        
        if strcmp(key, 'building')
            isBuilding = true;
            buildingType = value;
            break;
        end
    end
    
    if isBuilding
        buildingCount = buildingCount + 1;
        
        % 获取构成建筑物的节点
        nodeRefs = way.getElementsByTagName('nd');
        buildingCoords = [];
        
        for k = 0:nodeRefs.getLength()-1
            nodeRef = nodeRefs.item(k);
            nodeId = char(nodeRef.getAttribute('ref'));
            
            if isKey(nodeMap, nodeId)
                coord = nodeMap(nodeId);
                buildingCoords = [buildingCoords; coord];
            end
        end
        
        % 存储建筑物信息
        if ~isempty(buildingCoords)
            building.id = buildingCount;
            building.type = buildingType;
            building.latlon = buildingCoords;  % 保存原始经纬度 [lat, lon]
            
            buildings = [buildings; building];
        end
    end
end

fprintf('找到 %d 个建筑物\n', buildingCount);

if isempty(buildings)
    warning('未找到任何建筑物数据');
    return;
end

% 计算所有建筑物的边界，用于坐标转换
fprintf('正在计算坐标边界并转换为米制坐标...\n');
allLats = [];
allLons = [];
for i = 1:length(buildings)
    allLats = [allLats; buildings(i).latlon(:, 1)];
    allLons = [allLons; buildings(i).latlon(:, 2)];
end

minLat = min(allLats);
maxLat = max(allLats);
minLon = min(allLons);
maxLon = max(allLons);

% 使用地图中心点的纬度来计算经度转换系数
centerLat = (minLat + maxLat) / 2;
metersPerDegreeLat = 111320; % 1度纬度的米数
metersPerDegreeLon = 111320 * cos(deg2rad(centerLat)); % 1度经度的米数（考虑纬度影响）

fprintf('坐标转换参数:\n');
fprintf('  参考纬度: %.6f°\n', centerLat);
fprintf('  1度纬度 = %.2f 米\n', metersPerDegreeLat);
fprintf('  1度经度 = %.2f 米\n', metersPerDegreeLon);

% 转换所有建筑物坐标为米制坐标（相对于左下角原点）
for i = 1:length(buildings)
    latlon_coords = buildings(i).latlon;
    
    % 转换为米制坐标（相对于左下角原点）
    x_meters = (latlon_coords(:, 2) - minLon) * metersPerDegreeLon; % 经度差 -> 东向距离
    y_meters = (latlon_coords(:, 1) - minLat) * metersPerDegreeLat; % 纬度差 -> 北向距离
    
    % 存储米制坐标
    buildings(i).x_meters = x_meters;
    buildings(i).y_meters = y_meters;
    buildings(i).coordinates_meters = [x_meters, y_meters];
end

% 计算米制坐标范围
allXMeters = [];
allYMeters = [];
for i = 1:length(buildings)
    allXMeters = [allXMeters; buildings(i).x_meters];
    allYMeters = [allYMeters; buildings(i).y_meters];
end

minX = min(allXMeters);
maxX = max(allXMeters);
minY = min(allYMeters);
maxY = max(allYMeters);

% 保存数据到.mat文件
save(outputMatFile, 'buildings', 'minLat', 'minLon', 'maxLat', 'maxLon', ...
     'metersPerDegreeLat', 'metersPerDegreeLon', 'minX', 'maxX', 'minY', 'maxY');
fprintf('建筑物数据已保存到: %s\n', outputMatFile);

% 创建俯视图可视化（使用米制坐标）
figure('Name', 'OSM建筑物俯视图 (米制坐标)', 'Position', [100, 100, 1200, 800]);

% 设置地图边界（添加5%的边距）
xRange = maxX - minX;
yRange = maxY - minY;
xMargin = xRange * 0.05;
yMargin = yRange * 0.05;

xlim([minX - xMargin, maxX + xMargin]);
ylim([minY - yMargin, maxY + yMargin]);

hold on;

% 创建颜色映射
uniqueTypes = unique({buildings.type});
colors = hsv(length(uniqueTypes));
typeColorMap = containers.Map(uniqueTypes, num2cell(colors, 2));

% 绘制每个建筑物（使用米制坐标）
fprintf('正在绘制建筑物...\n');
for i = 1:length(buildings)
    building = buildings(i);
    
    if size(building.coordinates_meters, 1) >= 3
        % 确保多边形闭合
        coords_m = building.coordinates_meters;
        if ~isequal(coords_m(1,:), coords_m(end,:))
            coords_m = [coords_m; coords_m(1,:)];
        end
        
        % 获取建筑物类型对应的颜色
        if isKey(typeColorMap, building.type)
            color = typeColorMap(building.type);
        else
            color = [0.7, 0.7, 0.7]; % 默认灰色
        end
        
        % 绘制建筑物多边形（使用米制坐标）
        fill(coords_m(:, 1), coords_m(:, 2), color, ...
             'EdgeColor', 'black', 'LineWidth', 0.5, ...
             'FaceAlpha', 0.7);
    end
end

% 设置图形属性
grid on;
axis equal;
xlabel('东向距离 (米)');
ylabel('北向距离 (米)');
title(sprintf('OSM建筑物俯视图 (共%d个建筑) - 米制坐标', length(buildings)));

% 创建图例
if length(uniqueTypes) <= 10  % 只在类型不太多时显示图例
    legendHandles = [];
    legendLabels = {};
    
    for i = 1:length(uniqueTypes)
        type = uniqueTypes{i};
        color = typeColorMap(type);
        h = fill(NaN, NaN, color, 'EdgeColor', 'black', 'LineWidth', 0.5);
        legendHandles = [legendHandles, h];
        legendLabels{end+1} = sprintf('建筑类型: %s', type);
    end
    
    legend(legendHandles, legendLabels, 'Location', 'bestoutside');
end

hold off;

% 显示统计信息
fprintf('\n=== 建筑物统计信息 ===\n');
fprintf('总建筑物数量: %d\n', length(buildings));
fprintf('建筑物类型统计:\n');

typeCount = containers.Map();
for i = 1:length(buildings)
    type = buildings(i).type;
    if isKey(typeCount, type)
        typeCount(type) = typeCount(type) + 1;
    else
        typeCount(type) = 1;
    end
end

types = keys(typeCount);
for i = 1:length(types)
    fprintf('  %s: %d\n', types{i}, typeCount(types{i}));
end

fprintf('\n=== 坐标范围信息 ===\n');
fprintf('原始地图范围 (经纬度):\n');
fprintf('  纬度: %.6f° 到 %.6f°\n', minLat, maxLat);
fprintf('  经度: %.6f° 到 %.6f°\n', minLon, maxLon);

fprintf('转换后地图范围 (米制坐标):\n');
fprintf('  东向距离: %.2f 米 到 %.2f 米\n', minX, maxX);
fprintf('  北向距离: %.2f 米 到 %.2f 米\n', minY, maxY);
fprintf('  总宽度: %.2f 米\n', maxX - minX);
fprintf('  总高度: %.2f 米\n', maxY - minY);

% 可选：创建第二个图形显示建筑物分布密度（使用米制坐标）
figure('Name', '建筑物分布密度 (米制坐标)', 'Position', [200, 200, 800, 600]);

% 创建密度热图
nBins = 50;
[N, Xedges, Yedges] = histcounts2(allXMeters, allYMeters, nBins);

% 显示热图
imagesc(Xedges(1:end-1), Yedges(1:end-1), N');
colorbar;
xlabel('东向距离 (米)');
ylabel('北向距离 (米)');
title('建筑物分布密度热图 (米制坐标)');
axis xy;

fprintf('\n提取和可视化完成！\n');
fprintf('数据文件: %s\n', outputMatFile);

end

% 辅助函数：加载已保存的建筑物数据并重新绘制
function replotBuildingsFromMat(matFilePath)
% 从.mat文件重新绘制建筑物（米制坐标）
% 使用示例: replotBuildingsFromMat('buildings.mat');

if nargin < 1
    matFilePath = 'buildings.mat';
end

if ~exist(matFilePath, 'file')
    error('找不到文件: %s', matFilePath);
end

% 加载数据
load(matFilePath, 'buildings', 'minX', 'maxX', 'minY', 'maxY');
fprintf('从 %s 加载了 %d 个建筑物\n', matFilePath, length(buildings));

% 重新绘制（使用米制坐标）
figure('Name', '已保存的OSM建筑物俯视图 (米制坐标)');

% 设置地图边界
xRange = maxX - minX;
yRange = maxY - minY;
xMargin = xRange * 0.05;
yMargin = yRange * 0.05;

xlim([minX - xMargin, maxX + xMargin]);
ylim([minY - yMargin, maxY + yMargin]);

hold on;

% 创建颜色映射
uniqueTypes = unique({buildings.type});
colors = hsv(length(uniqueTypes));
typeColorMap = containers.Map(uniqueTypes, num2cell(colors, 2));

% 绘制每个建筑物
for i = 1:length(buildings)
    building = buildings(i);
    
    if size(building.coordinates_meters, 1) >= 3
        coords_m = building.coordinates_meters;
        if ~isequal(coords_m(1,:), coords_m(end,:))
            coords_m = [coords_m; coords_m(1,:)];
        end
        
        if isKey(typeColorMap, building.type)
            color = typeColorMap(building.type);
        else
            color = [0.7, 0.7, 0.7];
        end
        
        fill(coords_m(:, 1), coords_m(:, 2), color, ...
             'EdgeColor', 'black', 'LineWidth', 0.5, ...
             'FaceAlpha', 0.7);
    end
end

grid on;
axis equal;
xlabel('东向距离 (米)');
ylabel('北向距离 (米)');
title(sprintf('OSM建筑物俯视图 (共%d个建筑) - 米制坐标', length(buildings)));

hold off;

end

% 辅助函数：将经纬度坐标转换为米制坐标
function [x_meters, y_meters] = latlon2meters(lat, lon, ref_lat, ref_lon)
% 将经纬度转换为相对于参考点的米制坐标
% 输入:
%   lat, lon: 要转换的纬度和经度
%   ref_lat, ref_lon: 参考点的纬度和经度（通常是左下角）
% 输出:
%   x_meters: 东向距离（米）
%   y_meters: 北向距离（米）

metersPerDegreeLat = 111320;
metersPerDegreeLon = 111320 * cos(deg2rad(ref_lat));

x_meters = (lon - ref_lon) * metersPerDegreeLon;
y_meters = (lat - ref_lat) * metersPerDegreeLat;

end

% 辅助函数：将米制坐标转换回经纬度坐标
function [lat, lon] = meters2latlon(x_meters, y_meters, ref_lat, ref_lon)
% 将米制坐标转换回经纬度坐标
% 输入:
%   x_meters, y_meters: 米制坐标
%   ref_lat, ref_lon: 参考点的纬度和经度
% 输出:
%   lat, lon: 转换后的纬度和经度

metersPerDegreeLat = 111320;
metersPerDegreeLon = 111320 * cos(deg2rad(ref_lat));

lat = ref_lat + y_meters / metersPerDegreeLat;
lon = ref_lon + x_meters / metersPerDegreeLon;

end