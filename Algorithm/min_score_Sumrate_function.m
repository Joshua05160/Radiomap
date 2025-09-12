function score_R_sumrate = min_score_Sumrate_function(A,Power,pathloss, N_0)
   num_UAVs  = size(A,1);
   num_UGVs  = size(A,2);
   num_Slots = size(A,3);
   UGVrate = zeros(num_UGVs,1);
  % R_hat = zeros(num_UAVs,num_UGVs,num_Slots);
   for t = 1:num_Slots
       for n = 1:num_UGVs
           for m = 1:num_UAVs
               term1 = A(m,n,t)*Power(t,n)*pathloss(m,n,t);
               A_without_n = A(:,:,t);
               P_without_n =Power(t,:);
               pathloss_wihout_n = pathloss(m,:,t);
               pathloss_wihout_n(:,n) = [];
               P_without_n(:,n) =[];
               A_without_n(:,n) = [];
               Sum_a_pq = sum(A_without_n);
               term2 =sum(Sum_a_pq.*P_without_n.*pathloss_wihout_n);
               R_hat(m,n,t) = log2(term1+term2+N_0);
               R_tilde(m,n,t) = log2(term2+N_0);
               R_rate(m,n,t) =R_hat(m,n,t)-R_tilde(m,n,t);
           end
       end
   end
            for n = 1:num_UGVs
                UGVrate(n) = sum(sum(R_rate(:,n,:)));
            end

            score_R_sumrate = min(UGVrate);

end
               
               