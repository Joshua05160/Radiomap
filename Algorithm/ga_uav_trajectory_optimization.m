function [bestScore, scoreContainer, bestPosition, bestSumrate] = ga_uav_trajectory_optimization(M, T, BETA, A, Power, N_0, pathLossData, ALPHA, speedLimit)
    % 参数设置
    populationSize = 100;  % 种群大小
    numGenerations = 100;  % 迭代次数
    mutationRate = 0.1;  % 变异率
    crossoverRate = 0.8;  % 交叉率
    scoreContainer = zeros(numGenerations, 2);

    [population, velocities] = initialize_particles_multi(ALPHA(:,:,1), ALPHA(:,:,T), M, T, populationSize, ALPHA, speedLimit);
    % 初始化种群
    %population = initialize_population(ALPHA(:,:,1), ALPHA(:,:,T), M, T, populationSize, ALPHA, speedLimit);
    bestScore = -inf;
    bestPosition = population(:, :, :, 1);
    bestSumrate = 0;
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

    % 主循环
    for gen = 1:numGenerations
        scores = zeros(populationSize, 1);
        sumrates = zeros(populationSize, 1);

        % 计算适应度
        for i = 1:populationSize
            position = population(:, :, :, i);
            pathLoss = computePathloss(position, BETA, pathLossData);
            pathLoss = 10.^(-pathLoss / 10);
            scoreSumrate = score_Sumrate_function(A, Power, pathLoss, N_0);
            scoreSpeed = Speed_score_function(M, T, position, speedLimit);
            scoreAnticollision = collision_avoidance_fitness(position, building_heights_and_floors);
            angleScore = Angle_score_function(M, T, position, 50);

            scores(i) = 5 * scoreSpeed + 0.3 * scoreSumrate + 1.2 * scoreAnticollision + 5 * angleScore;
            sumrates(i) = scoreSumrate;
        end

        % 选择
        newPopulation = selection(population, scores, populationSize);

        % 交叉
        for i = 1:2:populationSize
            if rand < crossoverRate
                [child1, child2] = crossover(newPopulation(:, :, :, i), newPopulation(:, :, :, i+1));
                newPopulation(:, :, :, i) = child1;
                newPopulation(:, :, :, i+1) = child2;
            end
        end

        % 变异
        for i = 1:populationSize
            if rand < mutationRate
                newPopulation(:, :, :, i) = mutate(newPopulation(:, :, :, i), speedLimit);
            end
        end

        % 更新种群
        population = newPopulation;

        % 更新全局最佳位置
        [maxScore, maxIdx] = max(scores);
        if maxScore > bestScore
            bestScore = maxScore;
            bestPosition = population(:, :, :, maxIdx);
            bestSumrate = sumrates(maxIdx);
        end

        % 显示迭代信息
        fprintf('Generation %d: Best Score = %f, Best sumrate = %f\n', gen, bestScore,bestSumrate);
        scoreContainer(gen,:) = [bestScore bestSumrate];
    end
end

function population = initialize_population(startPosition, endPosition, M, T, populationSize, ALPHA, speedLimit)
    population = zeros(3, M, T, populationSize);
    for i = 1:populationSize
        population(:, :, 1, i) = startPosition;
        population(:, :, T, i) = endPosition;
        for t = 2:(T-1)
            for m = 1:M
                population(:, m, t, i) = startPosition + rand(3, 1) * speedLimit * (t - 1) / T;
            end
        end
    end
end

function newPopulation = selection(population, scores, populationSize)
    [~, sortedIdx] = sort(scores, 'descend');
    newPopulation = population(:, :, :, sortedIdx(1:populationSize));
end

function [child1, child2] = crossover(parent1, parent2)
    alpha = rand();
    child1 = alpha * parent1 + (1 - alpha) * parent2;
    child2 = (1 - alpha) * parent1 + alpha * parent2;
end

function mutatedPosition = mutate(position, speedLimit)
    [~, M, T] = size(position);
    for t = 2:(T-1)
        for m = 1:M
            mutation = (rand(3, 1) - 0.5) * 2 * speedLimit * 0.1;
            newPosition = position(:, m, t) + mutation;
            position(:, m, t) = newPosition;
        end
    end
    mutatedPosition = position;
end

% 请确保所有辅助函数（如 computePathloss、score_Sumrate_function、Speed_score_function、collision_avoidance_fitness、Angle_score_function 等）已实现，并与之前的代码保持一致。
