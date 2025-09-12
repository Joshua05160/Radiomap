function [Power_out, P_rate, P_sumrate]=Power_allocation_radiomap(M,N,T,A,gBestPosition,BETA,pathLossData,Power_l,N_0)
    cvx_clear
    cvx_begin quiet
  cvx_solver mosek
    variable P(T,N)

pathLoss = computePathloss(gBestPosition, BETA, pathLossData);
pathLoss = 10.^(-pathLoss/10);%pathLoss = computePathloss_withbuilding(positions(:,:,:,p), BETA, pathLossData,building_heights_and_floors);
R_hat = calculate_R_hat(A,P,pathLoss, N_0);
C_ub = calculate_C_ub_radiomap(A,P,Power_l,pathLoss,N_0);
    %psi_A = calculate_psi_A(A,A_l);
    A_state = sum(A,1);
    s=0;
for t=1:T
    for m=1:M
        for n=1:N
            s = s + R_hat(m,n,t)/log(2)-C_ub(m,n,t);
        end
    end
end
    maximize s
    
    subject to 
        for t=1:T
            for n =1:N
               0<=P(t,n)<=0.25;
            end
        end


        % for t=1:T
        %     for m=1:M
        %         for n=1:N
        %             0<=A(m,n,t)<=1;
        %             R_hat(m,n,t)/log(2)-C_ub(m,n,t)>=0.3*A(m,n,t);
        %         end
        %     end
        % end
        

    cvx_end
    Power_out = P;
    P_sumrate = cvx_optval;
    P_rate = R_hat/log(2)-C_ub;
    P_sumrate = sum(sum(sum(P_rate)));

end