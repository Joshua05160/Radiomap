function R_tilde_ub = calculate_R_tilde_ub(M,N,T,S,A,Distance_xy_l,Height,Power,L_0,N_0)
  % R_hat = zeros(num_UAVs,num_UGVs,num_Slots);
  
   for t = 1:T
       for n = 1:N
           for m = 1:M
               S_without_n = S(m,:,t);
               S_without_n(:,n) = [];
               
               
               H_without_n = ones(1,N-1);
               for ss = 1:N-1
                   H_without_n(ss) = Height(m,t);
               end
               A_without_n = A(:,:,t);
               P_without_n =Power(t,:);
               Distance_wihout_n = Distance_xy_l(m,:,t);
               Distance_wihout_n(:,n) = [];
               P_without_n(:,n) =[];
               A_without_n(:,n) = [];
               Sum_a_pq = sum(A_without_n);
               S_l =L_0*Sum_a_pq.*P_without_n./(Distance_wihout_n+H_without_n.^2);
               term1 = log2(sum(S_l)+N_0);
               term2 = sum(S_without_n-S_l)/(sum(S_l)+N_0);
               R_tilde_ub(m,n,t) = term1+term2/log(2);
           end
       end
   end