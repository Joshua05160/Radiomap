function score_R_sumrate= calculate_R_hat_R_tilde(A,Power,Distance,L_0,N_0)
   num_UAVs  = size(A,1);
   num_UGVs  = size(A,2);
   num_Slots = size(A,3);
  % R_hat = zeros(num_UAVs,num_UGVs,num_Slots);
   for t = 1:num_Slots
       for n = 1:num_UGVs
           for m = 1:num_UAVs
               term1 = A(m,n,t)*Power(t,n)*L_0/(Distance(m,n,t)^2);
               A_without_n = A(:,:,t);
               P_without_n =Power(t,:);
               Distance_wihout_n = Distance(m,:,t);
               Distance_wihout_n(:,n) = [];
               P_without_n(:,n) =[];
               A_without_n(:,n) = [];
               Sum_a_pq = sum(A_without_n);
               term2 =L_0*sum(Sum_a_pq.*P_without_n./(Distance_wihout_n.^2));
               R_hat(m,n,t) = log2(term1+term2+N_0);
               R_tilde(m,n,t) = log2(term2+N_0);
               R_rate(m,n,t) =R_hat(m,n,t)-R_tilde(m,n,t);

           end
       end
   end
      score_R_sumrate =sum(sum(sum(R_rate)));   
end       
               
               