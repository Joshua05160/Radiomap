function angle_score = Angle_score_function(M, T, positionsUAV, maxAngle)
    % 预分配空间
    penalty_sum = 0;

    % 计算每段轨迹之间的夹角
    for m = 1:M
        for t = 1:T-2
            % 计算向量
            vector1 = positionsUAV(:,m,t+1) - positionsUAV(:,m,t);
            vector2 = positionsUAV(:,m,t+2) - positionsUAV(:,m,t+1);
            
            % 计算向量的夹角
            dotProduct = dot(vector1, vector2);
            norm1 = norm(vector1);
            norm2 = norm(vector2);
            cosTheta = dotProduct / (norm1 * norm2);
            
            % 防止因浮点数运算误差导致cosTheta超出[-1,1]范围
            cosTheta = min(1, max(-1, cosTheta));
            
            % 计算角度（弧度制）
            angle = acos(cosTheta);
            
            % 将弧度转换为角度（度数制）
            angle_deg = rad2deg(angle);
            
            % 检查是否超过最大角度
            if angle_deg > maxAngle
                normalized_penalty = (angle_deg - maxAngle) / maxAngle;
                penalty_sum = penalty_sum + normalized_penalty;
            end
        end
    end
    
    % 适应度值 = 总惩罚值
    angle_score = -penalty_sum;
end
