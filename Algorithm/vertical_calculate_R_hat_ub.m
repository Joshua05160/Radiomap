function R_hat_ub = vertical_calculate_R_hat_ub(A,Power,Distance_xy,H_l,H,L_0,N_0)
   num_UAVs  = size(A,1);
   num_UGVs  = size(A,2);
   num_Slots = size(A,3);
  % R_hat = zeros(num_UAVs,num_UGVs,num_Slots);
   for t = 1:num_Slots
       for n = 1:num_UGVs
           for m = 1:num_UAVs
               term1 = A(m,n,t)*Power(t,n)*L_0/(Distance_xy(m,n,t)+H_l(m,t)^2);
               A_without_n = A(:,:,t);
               P_without_n =Power(t,:);
               Distance_xy_without_n = Distance_xy(m,:,t);
               Distance_xy_without_n(:,n) = [];
               %H_without_n = ones(1,num_UGVs-1);
               H_l_without_n = ones(1,num_UGVs-1);
               for ss = 1:num_UGVs-1
                   H_without_n(ss) = H(m,t);
                   H_l_without_n(ss) = H_l(m,t);
               end
               P_without_n(:,n) =[];
               A_without_n(:,n) = [];
               Sum_a_pq = sum(A_without_n);
               term2 =L_0*sum(Sum_a_pq.*P_without_n./(Distance_xy_without_n+H_l_without_n.^2));
               R_hat_ub_1(m,n,t) = log2(term1+term2+N_0);
               R_hat_ub_2(m,n,t) = A(m,n,t)*Power(t,n)*L_0/((Distance_xy(m,n,t)+H_l(m,t)^2))^2/(term1+term2+N_0)*(H(m,t)^2-H_l(m,t)^2);
               R_hat_ub_3(m,n,t) = L_0*sum((H_without_n.^2-H_l_without_n.^2).*Sum_a_pq.*P_without_n./(((Distance_xy_without_n+H_l_without_n.^2)).^2))/(term1+term2+N_0);
               R_hat_ub(m,n,t)=R_hat_ub_1(m,n,t)-R_hat_ub_2(m,n,t)/log(2)-R_hat_ub_3(m,n,t)/log(2);
           end
       end
   end
