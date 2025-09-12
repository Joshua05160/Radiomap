function [gBestScore,score_container,gBestPosition]=pso_uav_trajectory_optimization_noCM(M,T,BETA,A,Power,N_0,pathLossData,ALPHA,speedLimit,ini_method)
    % 参数设置
    numParticles = 100;  % 粒子数量
    numIterations =100;  % 迭代次数 % 速度限制 (m/s)
    score_container =zeros(numIterations,2);
    % 初始化粒子位置和速度
    if ini_method
        [positions, velocities] = initialize_particles(ALPHA(:,:,1), ALPHA(:,:,T), M, T, numParticles, ALPHA,speedLimit);
    else
        [positions, velocities] = initialize_particles_multi(ALPHA(:,:,1), ALPHA(:,:,T), M, T, numParticles, ALPHA, speedLimit);
    end
    % start_position 和 end_position 是 3x1 向量，表示 UAV 的起点和终点
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
    % 初始化粒子最佳位置和全局最佳位置
    pBestPositions = positions;
    pBestScores = -inf(numParticles, 1);
    gBestPosition = positions(:, :,:,1);
    gBestScore = -inf;

    % 主循环
    for iter = 1:numIterations
        for p = 1:numParticles
            % 更新速度
            velocities(:, :,:,p ) = update_velocity(velocities(:, :,:,p ), positions(:, :,:,p ), pBestPositions(:, :,:,p ), gBestPosition);
                % 固定起点和终点的速度为0
            velocities(:, :, 1, p ) = 0;
            velocities(:, :, T, p ) = 0;
            % 更新位置
            positions(:, :,:,p ) = positions(:, :,:,p) + velocities(:, :,:,p);
            pathLoss = computePathloss(positions(:,:,:,p), BETA, pathLossData);
            pathLoss = 10.^(-pathLoss/10);%pathLoss = computePathloss_withbuilding(positions(:,:,:,p), BETA, pathLossData,building_heights_and_floors);
            % % 计算适应度
            score_sumrate=score_Sumrate_function(A,Power,pathLoss, N_0);
            % 
            % heightweight = (2-0.7)*iter/numIterations+0.7;
            % [Distance,Distance_xy] = calculate_distance(gBestPosition_sca, squeeze(ALPHA(3, :,:)), BETA);
            % score_sumrate = calculate_R_hat_R_tilde(A,Power,Distance,10^(-3),N_0)
            score_speed = Speed_score_function(M,T,positions(:, :,:,p ), speedLimit);

            score_anticollision = collision_avoidance_fitness(positions(:, :,:,p), building_heights_and_floors);
            score_height_change = calculate_height_change(positions(:, :,:,p ));
            angle_score = Angle_score_function(M, T, positions(:, :,:,p ), 50);

            score = 5*score_speed+0.3*score_sumrate+1.2*score_anticollision+5*angle_score;
            % 1.2 3,1.2,1.1
            % 更新粒子最佳位置
            
            if score > pBestScores(p)
                pBestScores(p) = score;
                pBestPositions(:, :,:,p ) = positions(:, :,:,p );
            end
            
            % 更新全局最佳位置
            if score > gBestScore
                gBestScore = score;
                %score_container = [score,score_sumrate,score_speed,score_anticollision,score_height_change];
                gBestPosition = positions(:, :,:,p );
                gBestScore_sumrate = score_sumrate;
                gBestScore_speed = score_speed;
                gscore_anticollision = score_anticollision;
                gscore_height_change = score_height_change;
            end
        end

        % for p = 1:numParticles
        %     if rand < 0.1
        %         otherParticle = randi([1 numParticles]);
        %         positions(:, :, :, p) = crossover(positions(:, :, :, p), positions(:, :, :, otherParticle));
        %         velocities(:, :,:,p) = crossover(velocities(:, :,:,p), velocities(:, :, :, otherParticle));
        %     end
        % end
        % 
        % for p = 1:numParticles
        %     if rand < 0.1
        %         otherwaypoints = randi([1 T]);
        %         positions(:, :, otherwaypoints, p) = positions(:, :, otherwaypoints, p);
        %         velocities(:, :, otherwaypoints, p) = velocities(:, :, otherwaypoints, p);
        %     end
        % end
        % 
        % % 执行变异操作
        % for p = 1:numParticles
        %     if rand < 0.1
        %         positions(:, :, :, p) = mutate(positions(:, :, :, p), speedLimit);
        %         velocities(:, :, :, p) = mutate(velocities(:, :, :, p), speedLimit);
        % 
        %     end
        % end
        %显示迭代信息
        fprintf('Iteration %d: Best Score = %f, sumrate = %f, speed_score = %f, score_anticollision = %f, angle_score = %f\n', iter, gBestScore,gBestScore_sumrate,gBestScore_speed,gscore_anticollision,angle_score);
        score_container(iter,:) = [gBestScore gBestScore_sumrate];
    end
end

function new_velocity = update_velocity(velocity, position, pBestPosition, gBestPosition)
    % PSO 参数
    w = 0.5;  % 惯性权重
    c1 = 2;  % 个体学习因子
    c2 = 2;  % 群体学习因子
    % size(velocity)
    % size(position)
    % size(pBestPosition)
    % size(gBestPosition)
    % 更新速度
    inertia = w * velocity;
    cognitive = c1 * rand * (pBestPosition - position);
    social = c2 * rand * (gBestPosition - position);
    new_velocity = inertia + cognitive + social;
    % size(new_velocity)
end

function new_position = crossover(position1, position2)
    % 简单交叉操作
    alpha = rand();
    new_position = alpha * position1 + (1 - alpha) * position2;
end

function mutated_position = mutate(position, speedLimit)
    % 变异操作
    [~, M, T] = size(position);
    for t = 2:(T-1)
        for m = 1:M
            % 根据剩余步数计算每步最大移动距离
             % 随机变异
            mutation = (rand(3, 1) - 0.5) * 2 * speedLimit * 0.1; % 变异大小为10%的最大移动距离

            new_position = position(:, m, t) + mutation;
            
            position(:, m, t) = new_position;
        end
    end
    mutated_position = position;
end