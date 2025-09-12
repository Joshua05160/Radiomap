function B_ub = calculate_B_ub_radiomap(A,A_l,Power,pathLoss,N_0)
   num_UAVs  = size(A,1);
   num_UGVs  = size(A,2);
   num_Slots = size(A,3);
  % R_hat = zeros(num_UAVs,num_UGVs,num_Slots);
   for t = 1:num_Slots
       for n = 1:num_UGVs
           for m = 1:num_UAVs
               %term1 = A(m,n,t)*Power(t,n)*L_0/(Distance(m,n,t)^2);
               A_without_n = A(:,:,t);
               A_l_without_n = A_l(:,:,t);
               P_without_n =Power(t,:);
               pathLoss_wihout_n = pathLoss(m,:,t);
               pathLoss_wihout_n(:,n) = [];
               P_without_n(:,n) =[];
               A_without_n(:,n) = [];
               A_l_without_n(:,n) = [];
               Sum_a_pq = sum(A_without_n);
               Sum_a_l_pq = sum(A_l_without_n);
               summation_1 =sum(Sum_a_l_pq.*P_without_n.*pathLoss_wihout_n);
               term1 = log2(summation_1+N_0);
               summation_2 =sum((Sum_a_pq-Sum_a_l_pq).*P_without_n.*pathLoss_wihout_n);
               B_ub(m,n,t) = term1+(summation_2/((summation_1+N_0))*log(2));
           end
       end
   end