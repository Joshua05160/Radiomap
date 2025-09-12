function visualize_trajectory(UAV_trajecory, BETA,point_size,pathLossData, N_0)

load("building_data.mat");
% 提取建筑ID、X、Y和Z坐标
building_ids = unique(filtered_data.BuildingID);
num_buildings = length(building_ids);

% 创建图形窗口
figure;
hold on;
view(3);  % 设置为3D视图
grid on;  % 添加网格
axis equal;  % 坐标轴比例相同
% 绘制网格
grid_resolution = 5;
x_range = -100:grid_resolution:200;
y_range = -160:grid_resolution:160;
[X, Y] = meshgrid(x_range, y_range);
Z = zeros(size(X));
%mesh(X, Y, Z, 'FaceAlpha', 0, 'EdgeColor', 'k', 'LineStyle', '--')

sinr = show_sinrmap(BETA(:,:,5), 12,pathLossData,N_0);
sinr_map = sinr(:, :, 2);
sinr_map = flip(sinr_map,1);

 % imagesc(x_range, y_range, sinr_map);
 % set(gca, 'YDir', 'normal'); % 修正 Y 轴方向
 % colorbar; % 添加颜色条
 % caxis([0 6]);
 % colormap(jet);
   % 绘制无人机的轨迹和运动方向
% UAV_trajecory(1:2,:,:)=gBestPosition_sca;
% UAV_trajecory(3,:,:)=ALPHA(3,:,:);

ALPHA = UAV_trajecory;
    num_UAVs = size(ALPHA, 2);
    num_Vehicles = size(BETA, 2);
    num_points = size(BETA, 3);
    point_size = 10;
    T = num_points;
        for i = 1:num_UAVs
            % 从 ALPHA 矩阵中提取无人机的轨迹点
            trajectory = squeeze(ALPHA(:, i, :));
            x_values = trajectory(1, :);
            y_values = trajectory(2, :);
            z_values = trajectory(3, :);
            % 绘制轨迹（假设无人机在100米的高度上）
            plot3(x_values, y_values, z_values, 'LineWidth',2, 'Color', 'b');
            scatter3(x_values, y_values, z_values, point_size, 'filled', 'Marker', 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
            %start_end(i)=scatter3(ALPHA(1,i,1), ALPHA(2,i,1), 400, point_size, 'filled', 'Marker', 'hexagram', 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue','HandleVisibility', 'off');

            % 在起点处标记无人机编号
           text(x_values(1), y_values(1), z_values(1)+5, ['UAV ', num2str(i),' Start'], 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'blue');
           text(x_values(T), y_values(T), z_values(T)+5, ['UAV ', num2str(i),' End'], 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'blue');

            % 特殊标记起点
           plot3(x_values(1), y_values(1), z_values(1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
           plot3(x_values(T), y_values(T), z_values(T), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

            % 绘制无人机的运动方向
            %quiver3(x_values(1:end-1), y_values(1:end-1), ones(1, num_points-1) * 100, diff(x_values), diff(y_values), zeros(1, num_points-1), 'Color', 'red', 'LineWidth', 1.5);
        end
%     text(ALPHA(1,1)+30, ALPHA(2,1)-30, 400, ['UAV ', num2str(1),' (',num2str(round(ALPHA(1,1))),',',num2str(round(ALPHA(2,1))),',','100',')'], 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'blue');
%     text(ALPHA(1,2)+30, ALPHA(2,2)-30, 400, ['UAV ', num2str(2),' (',num2str(round(ALPHA(1,2))),',',num2str(round(ALPHA(2,2))),',','100',')'], 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'blue');
%     text(ALPHA(1,3)+30, ALPHA(2,3)-30, 400, ['UAV ', num2str(3),' (',num2str(round(ALPHA(1,3))),',',num2str(round(ALPHA(2,3))),',','100',')'], 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'blue');

    % 绘制车辆的轨迹和运动方向
    for i = 1:num_Vehicles
        % 从 BETA 矩阵中提取车辆的轨迹点
        trajectory = squeeze(BETA(:, i, :));
        x_values = trajectory(1, :);
        y_values = trajectory(2, :);
        
        % 绘制轨迹（假设车辆在地面上）
        plot3(x_values, y_values, zeros(1, num_points), 'LineWidth', 2, 'Color', 'red');
        scatter3(x_values, y_values, zeros(1, num_points), point_size, 'filled', 'Marker', 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
%         % 在起点处标记车辆编号
        text(x_values(1), y_values(1), 5, ['UGV ', num2str(i),' Start'], 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');
        text(x_values(T), y_values(T), 5, ['UGV ', num2str(i),' End'], 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'red');

        %quiver3(x_values(1)+30, y_values(1), 0, diff(x_values(1:2)), diff(y_values(1:2)), 0, 'Color', 'red', 'LineWidth', 3,'MaxHeadSize', 1);         
%         % 特殊标记起点
        plot3(x_values(1), y_values(1), 0, 'r','Marker','pentagram', 'MarkerSize', 8, 'MarkerFaceColor', 'magenta');
        plot3(x_values(T), y_values(T), 0, 'r', 'Marker','pentagram','MarkerSize', 8, 'MarkerFaceColor', 'magenta');

%         % 绘制车辆的运动方向
%         quiver3(x_values(20), y_values(20), 0, diff(x_values(1:2)), diff(y_values(1:2)), 0, 'Color', 'red', 'LineWidth', 1.5);
    end


% 定义颜色
top_face_color = hex2rgb('#C2C3BB'); 
side_face_color = hex2rgb('#F2F3F2');
side_line_color = hex2rgb('#9BA39A');

% 定义每个建筑的楼高和层数（根据建筑编号）
% 格式：building_heights_and_floors = [BuildingID, Height, Floors]
building_heights_and_floors = [
    1, 24, 8;
    2, 24, 8;
    3, 24, 8;
    4, 24, 8;
    5, 24, 8;
    6, 24, 8;
    7, 24, 8;
    8, 9, 3;
    9,36,12;
    10,42,14;
    11,48,16
];

% 遍历每个建筑并绘制其3D多边形
for i = 1:num_buildings
    building_data = filtered_data(filtered_data.BuildingID == building_ids(i), :);
    
    % 获取当前建筑的楼高和层数
    building_info = building_heights_and_floors(building_heights_and_floors(:, 1) == building_ids(i), :);
    building_height = building_info(2);
    floors = building_info(3);
    
    % 提取底部和顶部坐标
    X = [building_data.X; building_data.X(1)];
    Y = [building_data.Y; building_data.Y(1)];
    Z_bottom = zeros(size(X));
    Z_top = building_height * ones(size(X));
    
    % 绘制建筑的顶部
    fill3(X, Y, Z_top, side_face_color, 'FaceAlpha', 0.5);
    % 绘制建筑的底部
    fill3(X, Y, Z_bottom, side_face_color, 'FaceAlpha', 1);
    
    % 绘制每一层
    % for j = 1:floors
    %     z = (j - 1) * (building_height / floors);
    %     fill3(X, Y, z * ones(size(X)), side_face_color, 'FaceAlpha', 0, 'EdgeColor', side_line_color);
    % end
    
    % 绘制建筑的侧面
    for j = 1:length(X) - 1
        fill3([X(j), X(j+1), X(j+1), X(j)], [Y(j), Y(j+1), Y(j+1), Y(j)], ...
              [Z_bottom(j), Z_bottom(j+1), Z_top(j+1), Z_top(j)], side_face_color, 'FaceAlpha', 0.5);
    end
    
    % 计算建筑顶部的中心位置
    center_x = mean(X);
    center_y = mean(Y);
    center_z = max(Z_top);
    
    % 显示建筑编号
    % text(center_x, center_y, center_z, sprintf('%d', building_ids(i)), ...
    %      'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    %      'FontSize', 8, 'Color', 'blue');
end
    
% 圆柱1的参数
center_x = 31.063064409930682;
center_y = -15.243464178284768;
diameter1 = 3;
height1 = 2;
radius1 = diameter1 / 2;

% 圆柱2的参数
diameter2 = 2;
height2 = 1;
radius2 = diameter2 / 2;

% 圆台的参数
diameter3_bottom = 1;
diameter3_top = 0.5;
height3 = 6;
radius3_bottom = diameter3_bottom / 2;
radius3_top = diameter3_top / 2;

% 底部的高度
z_base_height = 48;

% 角度
theta = linspace(0, 2*pi, 8+1); % 8个顶点，正八边形

% 圆柱1的底部和顶部
x_cylinder1 = radius1 * cos(theta) + center_x;
y_cylinder1 = radius1 * sin(theta) + center_y;
z_cylinder1_bottom = z_base_height * ones(size(theta));
z_cylinder1_top = (z_base_height + height1) * ones(size(theta));

% 圆柱2的底部和顶部
x_cylinder2 = radius2 * cos(theta) + center_x;
y_cylinder2 = radius2 * sin(theta) + center_y;
z_cylinder2_bottom = (z_base_height + height1) * ones(size(theta));
z_cylinder2_top = (z_base_height + height1 + height2) * ones(size(theta));

% 圆台的底部和顶部
x_frustum_bottom = radius3_bottom * cos(theta) + center_x;
y_frustum_bottom = radius3_bottom * sin(theta) + center_y;
x_frustum_top = radius3_top * cos(theta) + center_x;
y_frustum_top = radius3_top * sin(theta) + center_y;
z_frustum_bottom = (z_base_height + height1 + height2) * ones(size(theta));
z_frustum_top = (z_base_height + height1 + height2 + height3) * ones(size(theta));

% 绘制圆柱和圆台
% 绘制圆柱1
fill3(x_cylinder1, y_cylinder1, z_cylinder1_bottom, top_face_color, 'FaceAlpha', 1)
fill3(x_cylinder1, y_cylinder1, z_cylinder1_top, top_face_color, 'FaceAlpha', 1)
for i = 1:length(theta)-1
    fill3([x_cylinder1(i), x_cylinder1(i+1), x_cylinder1(i+1), x_cylinder1(i)], ...
          [y_cylinder1(i), y_cylinder1(i+1), y_cylinder1(i+1), y_cylinder1(i)], ...
          [z_cylinder1_bottom(i), z_cylinder1_bottom(i+1), z_cylinder1_top(i+1), z_cylinder1_top(i)], side_face_color, 'FaceAlpha', 1)
end

% 绘制圆柱2
fill3(x_cylinder2, y_cylinder2, z_cylinder2_bottom, top_face_color, 'FaceAlpha', 1)
fill3(x_cylinder2, y_cylinder2, z_cylinder2_top, top_face_color, 'FaceAlpha', 1)
for i = 1:length(theta)-1
    fill3([x_cylinder2(i), x_cylinder2(i+1), x_cylinder2(i+1), x_cylinder2(i)], ...
          [y_cylinder2(i), y_cylinder2(i+1), y_cylinder2(i+1), y_cylinder2(i)], ...
          [z_cylinder2_bottom(i), z_cylinder2_bottom(i+1), z_cylinder2_top(i+1), z_cylinder2_top(i)], side_face_color, 'FaceAlpha', 1)
end

% 绘制圆台
fill3(x_frustum_bottom, y_frustum_bottom, z_frustum_bottom, top_face_color, 'FaceAlpha', 1)
fill3(x_frustum_top, y_frustum_top, z_frustum_top, top_face_color, 'FaceAlpha', 1)
for i = 1:length(theta)-1
    fill3([x_frustum_bottom(i), x_frustum_bottom(i+1), x_frustum_top(i+1), x_frustum_top(i)], ...
          [y_frustum_bottom(i), y_frustum_bottom(i+1), y_frustum_top(i+1), y_frustum_top(i)], ...
          [z_frustum_bottom(i), z_frustum_bottom(i+1), z_frustum_top(i+1), z_frustum_top(i)], side_face_color, 'FaceAlpha', 1)
end
% 设置图形属性
title('HITSZ', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (m)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('Height (m)', 'FontSize', 14, 'FontWeight', 'bold');

xlim([-70, 175]);
ylim([-130, 70]);

set(gca, 'LineWidth', 2); % 坐标轴加粗
set(gca, 'FontSize', 12); % 坐标轴上的数字变大
set(gca, 'FontWeight', 'bold'); % 坐标轴上的数字加粗

rotate3d on;

%     hg1   = gca;
%     hg2   = gca;
%     hg3   = gca;
%     hg4   = gca;
%     hg5   = gca;
%     hg6   = gca;
%     % 添加图例
%     %legend('UAV', 'Vehicle');
%     %% Design Different parts
%     %% Define design parameters
% D2R = pi/180;
% R2D = 180/pi;
% b   = 10;   % the length of total square cover by whole body of quadcopter in meter
% a   = b*2/3;   % the legth of small square base of quadcopter(b/4)
% H   = 1;  % hight of drone in Z direction (4cm)
% H_m = H+H/2; % hight of motor in z direction (5 cm)
% r_p = a*3/20;   % radius of propeller
% %% Conversions
% ro = 90*D2R;                   % angle by which rotate the base of quadcopter
% Ri = [cos(ro) -sin(ro) 0;
%       sin(ro) cos(ro)  0;
%        0       0       1];     % rotation matrix to rotate the coordinates of base 
% base_co = [-a/4  a/4 a/4 -a/4; % Coordinates of Base 
%            -a/2 -a/2 a/2 a/2;
%              0    0   0   0
%              ];
% base_co1 = [-a/4  -a/4 -a/4 -a/4; % Coordinates of Base 
%            -a/2 -a/2 a/2 a/2;
%              H    0   0   H
%              ];
% base_co2 = [-a/4  -a/4 a/4 a/4; % Coordinates of Base 
%            a/2 a/2 a/2 a/2;
%              0    H   H   0
%              ];
% base = Ri*base_co;             % rotate base Coordinates by 45 degree 
% base1 = Ri*base_co1; 
% base2 = Ri*base_co2; 
% to = linspace(0, 2*pi);
% xp = H/2*cos(to);
% yp = H/2*sin(to);
% zp = zeros(1,length(to));
% % design the base square
%  drone(1) = patch([base(1,:)],[base(2,:)],[base(3,:)],'r');
%  drone(2) = patch([base(1,:)],[base(2,:)],[base(3,:)+H],'r');
%  drone(3) = patch([base1(1,:)],[base1(2,:)],[base1(3,:)],'r');
%  drone(4) = patch([-base1(1,:)],[-base1(2,:)],[base1(3,:)],'r');
%  drone(5) = patch([base2(1,:)],[base2(2,:)],[base2(3,:)],'r');
%  drone(6) = patch([-base2(1,:)],[-base2(2,:)],[base2(3,:)],'r');
%  alpha(drone(1:6),0.7);
% % design 2 parpendiculer legs of quadcopter 
% [xcylinder ycylinder zcylinder] = cylinder([H/2 H/2]);
% %  drone(3) =  surface(b*zcylinder-b/2,ycylinder,xcylinder+H/2,'facecolor','b');
%  drone(7) =  surface(ycylinder+a/2-a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone(8) =  surface(ycylinder+a/2-a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone(9) =  surface(ycylinder-a/2+a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone(10) =  surface(ycylinder-a/2+a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  alpha(drone(7:10),0.6);
%  drone(11) =  surface((ycylinder/5+a/2-a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  drone(12) =  surface((ycylinder/5-a/2+a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  alpha(drone(11:12),0.6);
%  % % design 4 cylindrical motors 
% %  drone(5) = surface(xcylinder+b/2,ycylinder,H_m*zcylinder+H/2,'facecolor','r');
% %  drone(6) = surface(xcylinder-b/2,ycylinder,H_m*zcylinder+H/2,'facecolor','r');
% %  drone(7) = surface(xcylinder,ycylinder+b/2,H_m*zcylinder+H/2,'facecolor','r');
% %  drone(8) = surface(xcylinder,ycylinder-b/2,H_m*zcylinder+H/2,'facecolor','r');
% %  alpha(drone(5:8),0.7);
% % % design 4 propellers
%  drone(13)  = patch(xp+a/2-a*3/20,zp-a/4,yp+0.375 ,'p','LineWidth',0.3);
%  drone(14) = patch(xp-a/2+a*3/20,zp-a/4,yp+0.375,'p','LineWidth',0.3);
%  drone(15) = patch(xp+a/2-a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
%  drone(16) = patch(xp-a/2+a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
% 
%   drone(17)  = patch(xp+a/2-a*3/20,zp-a/4-a/15,yp+a/40 ,'p','LineWidth',0.3);
%  drone(18) = patch(xp-a/2+a*3/20,zp-a/4-a/15,yp+a/40,'p','LineWidth',0.3);
%  drone(19) = patch(xp+a/2-a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  drone(20) = patch(xp-a/2+a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  alpha(drone(13:20),1);
%  %drone1
%  drone1(1) = patch([base(1,:)],[base(2,:)],[base(3,:)],'r');
%  drone1(2) = patch([base(1,:)],[base(2,:)],[base(3,:)+H],'r');
%  drone1(3) = patch([base1(1,:)],[base1(2,:)],[base1(3,:)],'r');
%  drone1(4) = patch([-base1(1,:)],[-base1(2,:)],[base1(3,:)],'r');
%  drone1(5) = patch([base2(1,:)],[base2(2,:)],[base2(3,:)],'r');
%  drone1(6) = patch([-base2(1,:)],[-base2(2,:)],[base2(3,:)],'r');
%  alpha(drone1(1:6),0.7);
% % design 2 parpendiculer legs of quadcopter 
% [xcylinder ycylinder zcylinder] = cylinder([H/2 H/2]);
% %  drone(3) =  surface(b*zcylinder-b/2,ycylinder,xcylinder+H/2,'facecolor','b');
%  drone1(7) =  surface(ycylinder+a/2-a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone1(8) =  surface(ycylinder+a/2-a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone1(9) =  surface(ycylinder-a/2+a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone1(10) =  surface(ycylinder-a/2+a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  alpha(drone1(7:10),0.6);
%  drone1(11) =  surface((ycylinder/5+a/2-a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  drone1(12) =  surface((ycylinder/5-a/2+a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  alpha(drone1(11:12),0.6);
% 
%  drone1(13)  = patch(xp+a/2-a*3/20,zp-a/4,yp+0.375 ,'p','LineWidth',0.3);
%  drone1(14) = patch(xp-a/2+a*3/20,zp-a/4,yp+0.375,'p','LineWidth',0.3);
%  drone1(15) = patch(xp+a/2-a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
%  drone1(16) = patch(xp-a/2+a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
% 
%   drone1(17)  = patch(xp+a/2-a*3/20,zp-a/4-a/15,yp+a/40 ,'p','LineWidth',0.3);
%  drone1(18) = patch(xp-a/2+a*3/20,zp-a/4-a/15,yp+a/40,'p','LineWidth',0.3);
%  drone1(19) = patch(xp+a/2-a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  drone1(20) = patch(xp-a/2+a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  alpha(drone1(13:20),1);
% 
% %drone2
%  drone2(1) = patch([base(1,:)],[base(2,:)],[base(3,:)],'r');
%  drone2(2) = patch([base(1,:)],[base(2,:)],[base(3,:)+H],'r');
%  drone2(3) = patch([base1(1,:)],[base1(2,:)],[base1(3,:)],'r');
%  drone2(4) = patch([-base1(1,:)],[-base1(2,:)],[base1(3,:)],'r');
%  drone2(5) = patch([base2(1,:)],[base2(2,:)],[base2(3,:)],'r');
%  drone2(6) = patch([-base2(1,:)],[-base2(2,:)],[base2(3,:)],'r');
%  alpha(drone2(1:6),0.7);
% % design 2 parpendiculer legs of quadcopter 
% [xcylinder ycylinder zcylinder] = cylinder([H/2 H/2]);
% %  drone(3) =  surface(b*zcylinder-b/2,ycylinder,xcylinder+H/2,'facecolor','b');
%  drone2(7) =  surface(ycylinder+a/2-a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone2(8) =  surface(ycylinder+a/2-a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone2(9) =  surface(ycylinder-a/2+a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone2(10) =  surface(ycylinder-a/2+a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  alpha(drone2(7:10),0.6);
%  drone2(11) =  surface((ycylinder/5+a/2-a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  drone2(12) =  surface((ycylinder/5-a/2+a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  alpha(drone2(11:12),0.6);
% 
%  drone2(13)  = patch(xp+a/2-a*3/20,zp-a/4,yp+0.375 ,'p','LineWidth',0.3);
%  drone2(14) = patch(xp-a/2+a*3/20,zp-a/4,yp+0.375,'p','LineWidth',0.3);
%  drone2(15) = patch(xp+a/2-a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
%  drone2(16) = patch(xp-a/2+a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
% 
%   drone2(17)  = patch(xp+a/2-a*3/20,zp-a/4-a/15,yp+a/40 ,'p','LineWidth',0.3);
%  drone2(18) = patch(xp-a/2+a*3/20,zp-a/4-a/15,yp+a/40,'p','LineWidth',0.3);
%  drone2(19) = patch(xp+a/2-a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  drone2(20) = patch(xp-a/2+a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  alpha(drone2(13:20),1);
% 
%  %drone3
%   drone3(1) = patch([base(1,:)],[base(2,:)],[base(3,:)],'r');
%  drone3(2) = patch([base(1,:)],[base(2,:)],[base(3,:)+H],'r');
%  drone3(3) = patch([base1(1,:)],[base1(2,:)],[base1(3,:)],'r');
%  drone3(4) = patch([-base1(1,:)],[-base1(2,:)],[base1(3,:)],'r');
%  drone3(5) = patch([base2(1,:)],[base2(2,:)],[base2(3,:)],'r');
%  drone3(6) = patch([-base2(1,:)],[-base2(2,:)],[base2(3,:)],'r');
%  alpha(drone3(1:6),0.7);
% % design 2 parpendiculer legs of quadcopter 
% [xcylinder ycylinder zcylinder] = cylinder([H/2 H/2]);
% %  drone(3) =  surface(b*zcylinder-b/2,ycylinder,xcylinder+H/2,'facecolor','b');
%  drone3(7) =  surface(ycylinder+a/2-a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone3(8) =  surface(ycylinder+a/2-a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone3(9) =  surface(ycylinder-a/2+a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone3(10) =  surface(ycylinder-a/2+a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  alpha(drone3(7:10),0.6);
%  drone3(11) =  surface((ycylinder/5+a/2-a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  drone3(12) =  surface((ycylinder/5-a/2+a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  alpha(drone3(11:12),0.6);
% 
%  drone3(13)  = patch(xp+a/2-a*3/20,zp-a/4,yp+0.375 ,'p','LineWidth',0.3);
%  drone3(14) = patch(xp-a/2+a*3/20,zp-a/4,yp+0.375,'p','LineWidth',0.3);
%  drone3(15) = patch(xp+a/2-a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
%  drone3(16) = patch(xp-a/2+a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
% 
%   drone3(17)  = patch(xp+a/2-a*3/20,zp-a/4-a/15,yp+a/40 ,'p','LineWidth',0.3);
%  drone3(18) = patch(xp-a/2+a*3/20,zp-a/4-a/15,yp+a/40,'p','LineWidth',0.3);
%  drone3(19) = patch(xp+a/2-a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  drone3(20) = patch(xp-a/2+a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  alpha(drone3(13:20),1);
% 
% %% create a group object and parent surface
%   combinedobject1 = hgtransform('parent',hg1 );
%   combinedobject2 = hgtransform('parent',hg2 );
%   combinedobject3 = hgtransform('parent',hg3 );
%   combinedobject4 = hgtransform('parent',hg4 );
%   % set(drone,'parent',combinedobject1);
%   % set(drone1,'parent',combinedobject2);
%   % set(drone2,'parent',combinedobject3);
%   % set(drone3,'parent',combinedobject4);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%    
% 
%      % 
%      % translation1 = makehgtform('translate',...
%      %                           [BETA(1,1,1) BETA(2,1,1) 0]);
%      %  rotation3 = makehgtform('zrotate',pi/2+pi/4);
%      % set(combinedobject1,'matrix',...
%      %      translation1*rotation3);
%      % 
%      % 
%      % translation2 = makehgtform('translate',...
%      %                           [BETA(1,2,1) BETA(2,2,1) 0]);
%      %  rotation3 = makehgtform('zrotate',pi/2+pi/4);
%      % set(combinedobject2,'matrix',...
%      %      translation2*rotation3);
%      % 
%      % 
%      % translation3 = makehgtform('translate',...
%      %                           [BETA(1,3,1) BETA(2,3,1) 0]);
%      %  rotation3 = makehgtform('zrotate',pi/2+pi/4);
%      % set(combinedobject3,'matrix',...
%      %      translation3*rotation3);
% 
% 
%      % translation4 = makehgtform('translate',...
%      %                           [BETA(1,4,1) BETA(2,4,1) 0]);
%      %  rotation3 = makehgtform('zrotate',pi/2);
%      % set(combinedobject4,'matrix',...
%      %      translation4*rotation3);
%  drone5(1) = patch([base(1,:)],[base(2,:)],[base(3,:)],'r');
%  drone5(2) = patch([base(1,:)],[base(2,:)],[base(3,:)+H],'r');
%  drone5(3) = patch([base1(1,:)],[base1(2,:)],[base1(3,:)],'r');
%  drone5(4) = patch([-base1(1,:)],[-base1(2,:)],[base1(3,:)],'r');
%  drone5(5) = patch([base2(1,:)],[base2(2,:)],[base2(3,:)],'r');
%  drone5(6) = patch([-base2(1,:)],[-base2(2,:)],[base2(3,:)],'r');
%  alpha(drone5(1:6),0.7);
% % design 2 parpendiculer legs of quadcopter 
% [xcylinder ycylinder zcylinder] = cylinder([H/2 H/2]);
% %  drone(3) =  surface(b*zcylinder-b/2,ycylinder,xcylinder+H/2,'facecolor','b');
%  drone5(7) =  surface(ycylinder+a/2-a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone5(8) =  surface(ycylinder+a/2-a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone5(9) =  surface(ycylinder-a/2+a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone5(10) =  surface(ycylinder-a/2+a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  alpha(drone5(7:10),0.6);
%  drone5(11) =  surface((ycylinder/5+a/2-a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  drone5(12) =  surface((ycylinder/5-a/2+a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  alpha(drone5(11:12),0.6);
%  % % design 4 cylindrical motors 
% %  drone(5) = surface(xcylinder+b/2,ycylinder,H_m*zcylinder+H/2,'facecolor','r');
% %  drone(6) = surface(xcylinder-b/2,ycylinder,H_m*zcylinder+H/2,'facecolor','r');
% %  drone(7) = surface(xcylinder,ycylinder+b/2,H_m*zcylinder+H/2,'facecolor','r');
% %  drone(8) = surface(xcylinder,ycylinder-b/2,H_m*zcylinder+H/2,'facecolor','r');
% %  alpha(drone(5:8),0.7);
% % % design 4 propellers
%  drone5(13)  = patch(xp+a/2-a*3/20,zp-a/4,yp+0.375 ,'p','LineWidth',0.3);
%  drone5(14) = patch(xp-a/2+a*3/20,zp-a/4,yp+0.375,'p','LineWidth',0.3);
%  drone5(15) = patch(xp+a/2-a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
%  drone5(16) = patch(xp-a/2+a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
% 
%   drone5(17)  = patch(xp+a/2-a*3/20,zp-a/4-a/15,yp+a/40 ,'p','LineWidth',0.3);
%  drone5(18) = patch(xp-a/2+a*3/20,zp-a/4-a/15,yp+a/40,'p','LineWidth',0.3);
%  drone5(19) = patch(xp+a/2-a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  drone5(20) = patch(xp-a/2+a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  alpha(drone5(13:20),1);
% 
%    % combinedobject5 = hgtransform('parent',hg5 );
%    %   set(drone5,'parent',combinedobject5);
% 
%      %      translation5 = makehgtform('translate',...
%      %                           [BETA(1,5,1) BETA(2,5,1) 0]);
%      %  rotation3 = makehgtform('zrotate',pi/2);
%      % set(combinedobject5,'matrix',...
%      %      translation5*rotation3);
%        drone6(1) = patch([base(1,:)],[base(2,:)],[base(3,:)],'r');
%  drone6(2) = patch([base(1,:)],[base(2,:)],[base(3,:)+H],'r');
%  drone6(3) = patch([base1(1,:)],[base1(2,:)],[base1(3,:)],'r');
%  drone6(4) = patch([-base1(1,:)],[-base1(2,:)],[base1(3,:)],'r');
%  drone6(5) = patch([base2(1,:)],[base2(2,:)],[base2(3,:)],'r');
%  drone6(6) = patch([-base2(1,:)],[-base2(2,:)],[base2(3,:)],'r');
%  alpha(drone6(1:6),0.7);
% % design 2 parpendiculer legs of quadcopter 
% [xcylinder ycylinder zcylinder] = cylinder([H/2 H/2]);
% %  drone(3) =  surface(b*zcylinder-b/2,ycylinder,xcylinder+H/2,'facecolor','b');
%  drone6(7) =  surface(ycylinder+a/2-a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone6(8) =  surface(ycylinder+a/2-a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone6(9) =  surface(ycylinder-a/2+a*3/20,zcylinder-a/4-a/15,xcylinder+H/2-a/8,'facecolor','b') ; 
%  drone6(10) =  surface(ycylinder-a/2+a*3/20,zcylinder+a/4,xcylinder+H/2-a/8,'facecolor','b') ; 
%  alpha(drone6(7:10),0.6);
%  drone6(11) =  surface((ycylinder/5+a/2-a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  drone6(12) =  surface((ycylinder/5-a/2+a*3/20),a/2*zcylinder-a/4,xcylinder/5+H/2-a/8,'facecolor','k') ; 
%  alpha(drone6(11:12),0.6);
%  % % design 4 cylindrical motors 
% %  drone(5) = surface(xcylinder+b/2,ycylinder,H_m*zcylinder+H/2,'facecolor','r');
% %  drone(6) = surface(xcylinder-b/2,ycylinder,H_m*zcylinder+H/2,'facecolor','r');
% %  drone(7) = surface(xcylinder,ycylinder+b/2,H_m*zcylinder+H/2,'facecolor','r');
% %  drone(8) = surface(xcylinder,ycylinder-b/2,H_m*zcylinder+H/2,'facecolor','r');
% %  alpha(drone(5:8),0.7);
% % % design 4 propellers
%  drone6(13)  = patch(xp+a/2-a*3/20,zp-a/4,yp+0.375 ,'p','LineWidth',0.3);
%  drone6(14) = patch(xp-a/2+a*3/20,zp-a/4,yp+0.375,'p','LineWidth',0.3);
%  drone6(15) = patch(xp+a/2-a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
%  drone6(16) = patch(xp-a/2+a*3/20,zp+a/4,yp+0.375,'p','LineWidth',0.3);
% 
%   drone6(17)  = patch(xp+a/2-a*3/20,zp-a/4-a/15,yp+a/40 ,'p','LineWidth',0.3);
%  drone6(18) = patch(xp-a/2+a*3/20,zp-a/4-a/15,yp+a/40,'p','LineWidth',0.3);
%  drone6(19) = patch(xp+a/2-a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  drone6(20) = patch(xp-a/2+a*3/20,zp+a/4+a/15,yp+a/40,'p','LineWidth',0.3);
%  alpha(drone6(13:20),1);
% 
%    % combinedobject6 = hgtransform('parent',hg6 );
%    %   set(drone6,'parent',combinedobject6);
% 
%      %      translation6 = makehgtform('translate',...
%      %                           [BETA(1,6,1) BETA(2,6,1) 0]);
%      %  rotation3 = makehgtform('zrotate',pi/2);
%      % set(combinedobject6,'matrix',...
%      %      translation6*rotation3);
% 
% %     animation = drone_Animation(ALPHA(1,1,1), ALPHA(2,1,1),400,0,0,0);
% %     animation1 = drone_Animation(ALPHA(1,2,1), ALPHA(2,2,1),400,0,0,0);
% %     animation2 = drone_Animation(ALPHA(1,3,1), ALPHA(2,3,1),400,0,0,0);
% 
% 
%     %animation = drone_Animation(ALPHA(1,1,1), ALPHA(2,1,1),ALPHA(3,1,1),0,0,0);
%     %animation1 = drone_Animation(ALPHA(1,2,1), ALPHA(2,2,1),ALPHA(3,2,1),0,0,0);
%     % animation2 = drone_Animation(ALPHA(1,3,1), ALPHA(2,3,1),400,0,0,0);
hold off;
end

