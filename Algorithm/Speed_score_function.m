function speed_score = Speed_score_function(M,T,positionsUAV, speedLimit)
    % 计算每个时隙的距离
    for m=1:M
        for t=1:T-1
            distance(m,t)=pdist2(positionsUAV(:,m,t+1)',positionsUAV(:,m,t)','Euclidean');
        end
    end
        % 找到小于 0 的元素
    distance = 1.732*speedLimit-distance;
    elementsLessThanZero = distance(distance < 0);
    normalized_penalty = elementsLessThanZero / (1.732 * speedLimit);
    % 求和
    totalSum = sum(sum(normalized_penalty));
   
    % 适应度值 = 平均距离 + 惩罚项
    speed_score = totalSum;
end