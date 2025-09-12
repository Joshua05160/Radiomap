function R_hat = calculate_R_hat(A,Power,pathloss, N_0)
    num_UAVs  = size(A,1);
    num_UGVs  = size(A,2);
    num_Slots = size(A,3);
    %R_hat = zeros(num_UAVs, num_UGVs, num_Slots); % 初始化R_hat
    for t = 1:num_Slots
        for n = 1:num_UGVs
            for m = 1:num_UAVs
                term1 = A(m,n,t) * Power(t,n) * pathloss(m,n,t);
                A_without_n = A(:,:,t);
                P_without_n = Power(t,:);
                pathloss_without_n = pathloss(m,:,t);
                pathloss_without_n(:,n) = [];
                P_without_n(:,n) = [];
                A_without_n(:,n) = [];
                Sum_a_pq = sum(A_without_n);
                term2 = sum(Sum_a_pq .* P_without_n .* pathloss_without_n);
                R_hat(m,n,t) = log(term1 + term2 + N_0);
            end
        end
    end
   % R_hat = sum(sum(sum(R_hat))) / log(2) % 计算总和并转换单位
end
