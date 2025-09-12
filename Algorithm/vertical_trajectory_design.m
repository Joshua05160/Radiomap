function [H_out,verticle_Q_rate,verticle_Q_sumrate] =vertical_trajectory_design(M,N,T,A,Power,ALPHA_l,H_l,BETA,L_0,N_0,V_max,tau,d_min)
cvx_clear
cvx_begin quiet
% cvx_solver mosek
Z=2;
%ALPHA_initial= generate_line_trajectory(M, T, 600,350,10);
variable H(M,T)
variable B(M,N,T)
variable D(M,N,T)

[Distance,Distance_xy] = calculate_distance_2(ALPHA_l, H, BETA);
%[Distance_l,Distance_xy_l] = calculate_distance_2(ALPHA, H, BETA);

R_hat_ub = verticle_calculate_R_hat_ub(A,Power,Distance_xy,H_l,H,L_0,N_0);
R_tilde_ub = verticle_calculate_R_tilde_ub(M,N,T,B,A,Distance_xy,H_l,Power,L_0,N_0);

obj = 0;

     for t = 1:T
         for n = 1:N
             for m = 1:M
                  obj = obj+R_hat_ub(m,n,t)-R_tilde_ub(m,n,t);
             end
         end
     end
    maximize obj
    
    subject to
    
     for t = 1:T
         for n = 1:N
             for m = 1:M
                B(m,n,t)>=sum(A(:,n,t))*Power(t,n)*L_0/(norm(ALPHA_l(:,m,t)-BETA(:,n,t))^2+H_l(m,t)^2)-sum(A(:,n,t))*Power(t,n)*L_0/((norm(ALPHA_l(:,m,t)-BETA(:,n,t))^2+H_l(m,t)^2)^2)*(D(m,n,t)-norm(ALPHA_l(:,m,t)-BETA(:,n,t))^2-H_l(m,t)^2);
             end
         end
     end   
    H(:,1)==H_l(:,1);
    H(:,T)==H_l(:,T);
    %L(:,:,1)==ALPHA_initial(:,:,1);
%     L(:,:,1)==[200,500,800;100,100,100];
%     L(:,:,50)==[200,500,800;100,100,100];
     for t = 1:T
         for n = 1:N
             for m = 1:M
                 D(m,n,t)<= (ALPHA_l(1,m,t)-BETA(1,n,t))^2+2*H_l(m,t)*(H(m,t)-H_l(m,t));
                 %R_hat_ub(m,n,t)-R_tilde_ub(m,n,t)>=0.2*A(m,n,t);
             end
         end
     end
     for t = 1:T
             for m = 1:M
                 H(m,t)>=10;
                 %R_hat_ub(m,n,t)-R_tilde_ub(m,n,t)>=0.2*A(m,n,t);
             end

     end
%     for n=1:N
%          average_sumrate = 0;
%          for t=1:T
%              for m=1:M
%                  average_sumrate = average_sumrate+R_hat_ub(m,n,t)-R_tilde_ub(m,n,t);
%              end
%          end
%          average_sumrate/T >=0.1;
%     end
             
     
%      for t=1:T
%          for m = 1:M 
%             ALPHA_without_m = L(:,:,t);
%             ALPHA_l_without_m = ALPHA_l(:,:,t);
%             ALPHA_without_m (:,m)= [];
%             ALPHA_l_without_m (:,m)= [];
%             for p = 1:M-1
%                 d_min^2 <= norm(ALPHA_l(:,m,t)-ALPHA_l_without_m(:,p))^2+2*((ALPHA_l(:,m,t)-ALPHA_l_without_m(:,p))')*(L(:,m,t)-ALPHA_without_m(:,p)-(ALPHA_l(:,m,t)-ALPHA_l_without_m(:,p)));
%             end
%          end
%      end
          
     for t=1:T-1
         for m =1:M
         (ALPHA_l(1,m,t+1)-ALPHA_l(1,m,t))^2+(ALPHA_l(2,m,t+1)-ALPHA_l(2,m,t))^2+(H(m,t+1)-H(m,t))^2 <= V_max*tau;
         end
     end
     
cvx_end
     
     H_out = H;
     verticle_Q_rate = R_hat_ub-R_tilde_ub;
     verticle_Q_sumrate = obj;
end