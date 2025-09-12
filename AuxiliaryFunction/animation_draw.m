function animation_draw(ALPHA, BETA,Communication_A,point_size,file_path)
    num_UAVs = size(ALPHA, 2);
    num_Vehicles = size(BETA, 2);
    num_points = size(BETA, 3);
    %angle = pi/100;

    load("building_data.mat");
% 提取建筑ID、X、Y和Z坐标
    building_ids = unique(filtered_data.BuildingID);
    num_buildings = length(building_ids);
    % 创建新图形
    h=figure;
    hold on;
    view(3);  % 设置为3D视图
    grid on;  % 添加网格
    axis equal;  % 坐标轴比例相同
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	% 关闭相关的警告提示（因为调用了非公开接口）
    jFrame = get(h,'JavaFrame');	% 获取底层 Java 结构相关句柄吧
    pause(0.1);					% 在 Win 10，Matlab 2017b 环境下不加停顿会报 Java 底层错误。各人根据需要可以进行实验验证
    set(jFrame,'Maximized',1);	%设置其最大化为真（0 为假）
    pause(0.1);					% 个人实践中发现如果不停顿，窗口可能来不及变化，所获取的窗口大小还是原来的尺寸。各人根据需要可以进行实验验证
    warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		% 打开相关警告设置

 %ALPHA = UAV_trajecory;
    num_UAVs = size(ALPHA, 2);
    num_Vehicles = size(BETA, 2);
    num_points = size(BETA, 3);
    point_size = 10;
    UAVplot = zeros(num_UAVs,1);
    UGVplot = zeros(num_Vehicles,1);
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
           UAVplot(i) = plot3(x_values(1), y_values(1), z_values(1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
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
       UGVplot(i)=plot3(x_values(1), y_values(1), 0, 'r','Marker','pentagram', 'MarkerSize', 8, 'MarkerFaceColor', 'magenta');
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

hold on;


interferenceplot=zeros(num_UAVs,1);
communicationplot=zeros(num_UAVs,1);

    for m=1:num_UAVs
        for n=1:num_Vehicles
            if Communication_A(m,n,1)==1              
              communicationplot(m)=quiver3(BETA(1,n,1),BETA(2,n,1),0,ALPHA(1,m,1)-BETA(1,n,1),ALPHA(2,m,1)-BETA(2,n,1),ALPHA(3,m,1),0,'k-','LineWidth', 2);
              interferedUAV = [1,2];
              interferedUAV(m)=[];
              
              interferenceplot(m)=plot3([BETA(1,n,1),ALPHA(1,interferedUAV,1)],[BETA(2,n,1), ALPHA(2,interferedUAV,1)],[0,ALPHA(3,interferedUAV,1)],'k--','LineWidth', 1);
              %interferenceplot(m,2)=plot3([BETA(1,n,1),ALPHA(1,interferedUAV(2),1)],[BETA(2,n,1), ALPHA(2,interferedUAV(2),1)],[0,400],'k--','LineWidth', 1);
              %h_1plot(j)=plot3([x_path(1, i),uav(j, 1)],[ y_path(1, i), uav(j, 2)],[0,2*uav(j,3)],'k-'); 
            end
        end
    end
%legend([communicationplot(1),interferenceplot(1)],'Communication link','Interference link','AutoUpdate','off');
%legend('Communication link','Interference link','AutoUpdate','off');

new_ALPHA = interpolate_trajectory(ALPHA, 200);
new_BETA = interpolate_trajectory(BETA, 200);
% new_ALPHA=ALPHA;
% new_BETA = BETA;

pause(3);
frame = getframe(gcf);
imind = frame2im(frame);
[imind, cm] = rgb2ind(imind,256);
imwrite(imind,cm, file_path, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
for i = 2:200
    % delete(animation);
    % delete(animation1);
    % delete(animation2);
    for m=1:num_UAVs
        set(UAVplot(m),'XData', new_ALPHA(1,m,i), 'YData',new_ALPHA(2,m,i),'ZData',new_ALPHA(3,m,i));
    end
    for n=1:num_Vehicles
        set(UGVplot(n),'XData', new_BETA(1,n,i), 'YData',new_BETA(2,n,i),'ZData',0);
    end
    % animation = drone_Animation(new_ALPHA(1,1,i), new_ALPHA(2,1,i),400,0,0,0);
    % animation1 = drone_Animation(new_ALPHA(1,2,i), new_ALPHA(2,2,i),400,0,0,0);
    % animation2 = drone_Animation(new_ALPHA(1,3,i), new_ALPHA(2,3,i),400,0,0,0);

      
        dex=ceil(i/10);
%         dex=i;
       %title( ['Time slot=',num2str(dex/2+1)]);
       %text(300,750,['Time slot= ',num2str(dex)]);
    for m=1:num_UAVs
        for n=1:num_Vehicles
            if Communication_A(m,n,dex)==1              
              interferedUAV = [1,2];
              interferedUAV(m)=[];
                 set(interferenceplot(m),'XData', [new_BETA(1,n,i),new_ALPHA(1,interferedUAV,i)], 'YData',[new_BETA(2,n,i), new_ALPHA(2,interferedUAV,i)],'ZData',[0,new_ALPHA(3,interferedUAV,i)]);
                 %set(interferenceplot(m,2),'XData', [new_BETA(1,n,i),new_ALPHA(1,interferedUAV(2),i)], 'YData',[new_BETA(2,n,i), new_ALPHA(2,interferedUAV(2),i)],'ZData',[0,400]);
                 set(communicationplot(m),'Xdata',new_BETA(1,n,i),'Ydata',new_BETA(2,n,i),'Zdata',0,'Udata',new_ALPHA(1,m,i)-new_BETA(1,n,i),'Vdata',new_ALPHA(2,m,i)-new_BETA(2,n,i),'Wdata',new_ALPHA(3,m,i));
              %h_1plot(j)=plot3([x_path(1, i),uav(j, 1)],[ y_path(1, i), uav(j, 2)],[0,2*uav(j,3)],'k-'); 
            end
        end
    end

        pause(0.1);
frame = getframe(gcf);
imind = frame2im(frame);
[imind, cm] = rgb2ind(imind,256);
 imwrite(imind, cm, file_path, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
% legend([UGV_paths(1),interferenceplot(1),communicationplot(1)],'UGV path','Interference link','Communication link','Location', 'Northwest','AutoUpdate','off');
hold on;

end
hold off
end

function rgb = hex2rgb(hex)
    % Convert a hex color string to an RGB vector
    hex = hex(2:end);  % Remove the '#'
    r = hex2dec(hex(1:2))/255;
    g = hex2dec(hex(3:4))/255;
    b = hex2dec(hex(5:6))/255;
    rgb = [r, g, b];
end
