function [Sche, Optimal_val,R_out,psi_A_out]=communication_scheduling_radiomap(M,N,T,A_l,Power,gBestPosition,BETA,pathLossData,N_0,state)
    cvx_clear
    cvx_begin quiet
    cvx_solver mosek
    variable A(M,N,T)
    variable fairness
    
    pathLoss = computePathloss(gBestPosition, BETA, pathLossData);
    pathLoss = 10.^(-pathLoss/10);
    R_hat = calculate_R_hat(A,Power,pathLoss, N_0);
    B_ub = calculate_B_ub_radiomap(A,A_l,Power,pathLoss,N_0);
    psi_A = A.*A_l.*2-A-A_l.^2;
    s=0;
     for t = 1:T
         for n = 1:N
             for m = 1:M
                  s=s+R_hat(m,n,t)/log(2)-B_ub(m,n,t)+0.3*psi_A(m,n,t);
                 
             end
         end
     end
    maximize fairness

    subject to 
        for t=1:T
            for n =1:N
                0<=sum(A(:,n,t))<=1;
            end
        end

        if state
            for t=1:T
                    A(1,3,t)==1;
            end
        end
        for t=1:T
            for m=1:M
                0<=sum(A(m,:,t))<=1;
            end
        end

        for t=1:T
            for m=1:M
                for n=1:N
                    0<=A(m,n,t)<=1;
                   % R_hat(m,n,t)/log(2)-B_ub(m,n,t)>=0.3*A(m,n,t);
                end
            end
        end
        
    for n=1:N
         average_sumrate = 0;
         for t=1:T
             for m=1:M
                 average_sumrate = average_sumrate+R_hat(m,n,t)/log(2)-B_ub(m,n,t);
             end
         end
         average_sumrate/T >=fairness;
    end

    cvx_end
    Sche = A;
    Optimal_val = cvx_optval;
    R_out = R_hat/log(2)-B_ub;
    psi_A_out =  sum(sum(sum(psi_A)));
end