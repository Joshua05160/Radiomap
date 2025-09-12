T = 20; % time slots
M = 2; % UAV number
N = 4; % UGV number
A_l = 0.1*ones(M,N,T);
d_min = 10; %anti-collision distance /m
N_0 = 10^(-120/10);%noise power /dB
L_0 = 10^(-30/10);

BETA = zeros(2,N,T);
ALPHA = zeros(3,M,T);
Scheduling = zeros(M,N,T);
Power=0.5*ones(T,N);


A_scapsocm = zeros(M,N,T,100);
A_psocm = zeros(M,N,T,100);
A_pso = zeros(M,N,T,100);
A_ga = zeros(M,N,T,100);
A_close = zeros(M,N,T,100);
A_halfpower = zeros(M,N,T,100);

P_scapsocm = zeros(T,N,100);
P_psocm = zeros(T,N,100);
P_pso = zeros(T,N,100);
P_ga = zeros(T,N,100);
P_close = zeros(T,N,100);
P_halfpower = zeros(T,N,100);

%%%SCA-PSO-CM Solve
gBestPosition = ALPHA;
increment_total =100;
Power_out = 0.5*ones(T,N);
Sum_ini  = 0;
i = 0;
H_out = squeeze(ALPHA(3,:,:));
%gBestPosition = pisition_container_scapsocm(:,:,:,10);
position_container_scapsocm= zeros(3,M,T,100);
score_container_scapsocm = zeros(100,2,100);
while i<5
    i = i+1 ;
    %Solve Scheduling Problem
    increment = 100;
    Optimal_val_ini = -100;
    A_l = 0.1*ones(M,N,T);
    while increment>0.01
    [Sche, Optimal_val,R_out,psi_A_out]=communication_scheduling_radiomap(M,N,T,A_l,Power_out,gBestPosition,BETA,pathLossData,N_0,0);
    A_l = Sche;
    increment = Optimal_val - Optimal_val_ini;
    Optimal_val_ini = Optimal_val
    end
    Optimal_val
    A_scapsocm(:,:,:,i) = Sche;
    % Power_l = 0.1*ones(T,N);
    % power_container = zeros(T,N,20);
    % P_sumrate_ini = 0;
    % increment = 100;
    % j =1;
    % while increment>0.1
    % [Power_out, P_rate, P_sumrate]=Power_allocation_radiomap(M,N,T,round(Sche),gBestPosition,BETA,pathLossData,Power_l,N_0);
    % power_container(:,:,j)=Power_out;
    % j=j+1;
    % Power_l = Power_out;
    % increment = P_sumrate - P_sumrate_ini;
    % P_sumrate_ini = P_sumrate
    % end
    % P_sumrate
    % P_scapsocm(:,:,i) = Power_out;
    %迭代求解trajectory design
    for xx = i:5
        gBestPosition_sca = ALPHA(1:2,:,:);
        Q_sumrate_ini = -100;
        increment = 100;
        while increment>0.1
        [iter_position_sca,Q_rate,Q_sumrate] = trajectory_design(M,N,T,round(Sche),Power_out,gBestPosition_sca,H_out,BETA,L_0,N_0,20,20,d_min);
        gBestPosition_sca = iter_position_sca;
        increment = Q_sumrate - Q_sumrate_ini;
        Q_sumrate_ini = Q_sumrate
        end
        Q_sumrate
    
        %Sovle Height Problem
        H_l = squeeze(ALPHA(3,:,:));
        increment = 100;
        verticl_ini = -100;
        while increment>0.1
        [H_out,verticle_Q_rate,verticle_Q_sumrate] =vertical_trajectory_design(M,N,T,round(Sche),Power_out,gBestPosition_sca,H_l,BETA,L_0,N_0,20,20,d_min);
        H_l = H_out;
        increment = verticle_Q_sumrate - verticl_ini;
        verticl_ini = verticle_Q_sumrate
        end
    
        gBestPosition(1:2,:,:) = gBestPosition_sca;
        gBestPosition(3,:,:) = H_out;
    end


    % gBestPosition(1:2,:,:) = gBestPosition_sca;
    % gBestPosition(3,:,:) = H_out;
    % for j =1:5
    %     [gBestScore,container,iter_position]=pso_uav_trajectory_optimization(M,T,BETA,round(Sche),Power_out,N_0,pathLossData,gBestPosition,20/sqrt(3),1);
    %     position_container_scapsocm(:,:,:,5*(i-1)+j)= iter_position;
    %     score_container_scapsocm(:,:,5*(i-1)+j) = container;
    %     gBestPosition = iter_position;
    %     gBestPosition = smooth_trajectory(gBestPosition, 3);
    %     fprintf("这是SCA-PSO-CM第%d大循环 第%d小循环",i,j);
    % end

end



gBestPosition = ALPHA;
increment_total =100;
Power_out = 0.5*ones(T,N);
Sum_ini  = 0;
i = 0;
H_out = squeeze(ALPHA(3,:,:));
%gBestPosition = pisition_container_scapsocm(:,:,:,10);
position_container_scapsocm= zeros(3,M,T,100);
score_container_scapsocm = zeros(100,2,100);
while i<5
    i=i+1;
    Power_l = 0.1*ones(T,N);
    power_container = zeros(T,N,20);
    P_sumrate_ini = 0;
    increment = 100;
    j =1;
    while increment>0.1
    [Power_out, P_rate, P_sumrate]=Power_allocation_radiomap(M,N,T,round(Sche),gBestPosition,BETA,pathLossData,Power_l,N_0);
    power_container(:,:,j)=Power_out;
    j=j+1;
    Power_l = Power_out;
    increment = P_sumrate - P_sumrate_ini;
    P_sumrate_ini = P_sumrate
    end
    P_sumrate
    P_scapsocm(:,:,i) = Power_out;

    gBestPosition(1:2,:,:) = gBestPosition_sca;
    gBestPosition(3,:,:) = H_out;
    for j =1:5
        [gBestScore,container,iter_position]=pso_uav_trajectory_optimization(M,T,BETA,round(Sche),Power_out,N_0,pathLossData,gBestPosition,20/sqrt(3),1);
        position_container_scapsocm(:,:,:,5*(i-1)+j)= iter_position;
        score_container_scapsocm(:,:,5*(i-1)+j) = container;
        gBestPosition = iter_position;
        gBestPosition = smooth_trajectory(gBestPosition, 3);
        fprintf("这是SCA-PSO-CM第%d大循环 第%d小循环",i,j);
    end

end


%%% PSO-CM Solve
iter_position_pso_cm = ALPHA;
increment_total =100;
Power_out = 0.5*ones(T,N);
Sum_ini  = 0;
i = 0;
H_out = squeeze(ALPHA(3,:,:));
%gBestPosition = pisition_container_scapsocm(:,:,:,10);
pisition_container_pso_cm= zeros(3,M,T,100);
score_container_pso_cm = zeros(100,2,100);
while i<5
    i = i+1 ;
    %Solve Scheduling Problem
    increment = 100;
    Optimal_val_ini = -100;
    A_l = 0.1*ones(M,N,T);
    while increment>0.01
    [Sche, Optimal_val,R_out,psi_A_out]=communication_scheduling_radiomap(M,N,T,A_l,Power_out,iter_position_pso_cm,BETA,pathLossData,N_0,0);
    A_l = Sche;
    increment = Optimal_val - Optimal_val_ini;
    Optimal_val_ini = Optimal_val;
    end
    Optimal_val

    % Power_l = 0.1*ones(T,N);
    % power_container = zeros(T,N,20);
    % P_sumrate_ini = 0;
    % increment = 100;
    % j =1;
    % while increment>0.1
    % [Power_out, P_rate, P_sumrate]=Power_allocation_radiomap(M,N,T,round(Sche),iter_position_pso_cm,BETA,pathLossData,Power_l,N_0);
    % power_container(:,:,j)=Power_out;
    % j=j+1;
    % Power_l = Power_out;
    % increment = P_sumrate - P_sumrate_ini;
    % P_sumrate_ini = P_sumrate
    % end
    % P_sumrate
    % P_psocm(:,:,i) = Power_out;

    A_psocm(:,:,:,i) = Sche;
    gBestPosition_pso_cm=ALPHA;
    for j =1:5
        [gBestScore_pso_cm,container_pso_cm,iter_position_pso_cm]=pso_uav_trajectory_optimization(M,T,BETA,round(Sche),Power_out,N_0,pathLossData,gBestPosition_pso_cm,20/sqrt(3),0);
        pisition_container_pso_cm(:,:,:,5*(i-1)+j)= iter_position_pso_cm;
        score_container_pso_cm(:,:,5*(i-1)+j) = container_pso_cm;
        % gBestPosition_pso = iter_position_pso;
        % 
        % gBestPosition_pso = smooth_trajectory(gBestPosition_pso, 3);
        fprintf("这是PSO-CM第%d大循环 第%d小循环",i,j);

    end

end



%%% PSO Solve
gBestPosition_pso = ALPHA;
increment_total =100;
Power_out = 0.5*ones(T,N);
Sum_ini  = 0;
i = 0;
%gBestPosition = pisition_container_scapsocm(:,:,:,10);
pisition_container_pso= zeros(3,M,T,100);
score_container_pso = zeros(100,2,100);
while i<5
    i = i+1 ;
    %Solve Scheduling Problem
    increment = 100;
    Optimal_val_ini = -100;
    A_l = 0.1*ones(M,N,T);
    while increment>0.01
    [Sche, Optimal_val,R_out,psi_A_out]=communication_scheduling_radiomap(M,N,T,A_l,Power_out,gBestPosition_pso,BETA,pathLossData,N_0,0);
    A_l = Sche;
    increment = Optimal_val - Optimal_val_ini;
    Optimal_val_ini = Optimal_val;
    end
    Optimal_val
    A_pso(:,:,:,i) = Sche;

    % Power_l = 0.1*ones(T,N);
    % power_container = zeros(T,N,20);
    % P_sumrate_ini = 0;
    % increment = 100;
    % j =1;
    % while increment>0.1
    % [Power_out, P_rate, P_sumrate]=Power_allocation_radiomap(M,N,T,round(Sche),gBestPosition_pso,BETA,pathLossData,Power_l,N_0);
    % power_container(:,:,j)=Power_out;
    % j=j+1;
    % Power_l = Power_out;
    % increment = P_sumrate - P_sumrate_ini;
    % P_sumrate_ini = P_sumrate
    % end
    % P_sumrate
    %     P_pso(:,:,i) = Power_out;

    gBestPosition_pso=ALPHA;
    pisition_container_pso = zeros(3,M,T,20);
    score_container_pso = zeros(100,2,20);
    for j =1:5
        [gBestScore_pso,container_pso,iter_position_pso]=pso_uav_trajectory_optimization_noCM(M,T,BETA,round(Sche),Power_out,N_0,pathLossData,gBestPosition_pso,20/sqrt(3),0);
        pisition_container_pso(:,:,:,5*(i-1)+j)= iter_position_pso;
        score_container_pso(:,:,5*(i-1)+j) = container_pso;
        fprintf("这是PSO第%d大循环 第%d小循环",i,j);

    end

    gBestPosition_pso = iter_position_pso;
end
%%%





%%% GM Solve
gBestPosition_ga = ALPHA;
increment_total =100;
Power_out = 0.5*ones(T,N);
Sum_ini  = 0;
i = 0;
%gBestPosition = pisition_container_scapsocm(:,:,:,10);
pisition_container_ga= zeros(3,M,T,100);
score_container_ga = zeros(100,2,100);
while i<5
    i = i+1 ;
    %Solve Scheduling Problem
    increment = 100;
    Optimal_val_ini = -100;
    A_l = 0.1*ones(M,N,T);
    while increment>0.01
    [Sche, Optimal_val,R_out,psi_A_out]=communication_scheduling_radiomap(M,N,T,A_l,Power_out,gBestPosition_ga,BETA,pathLossData,N_0,0);
    A_l = Sche;
    increment = Optimal_val - Optimal_val_ini;
    Optimal_val_ini = Optimal_val;
    end
    Optimal_val
    A_ga(:,:,:,i) = Sche;

    % Power_l = 0.1*ones(T,N);
    % power_container = zeros(T,N,20);
    % P_sumrate_ini = 0;
    % increment = 100;
    % j =1;
    % while increment>0.1
    % [Power_out, P_rate, P_sumrate]=Power_allocation_radiomap(M,N,T,round(Sche),gBestPosition_ga,BETA,pathLossData,Power_l,N_0);
    % power_container(:,:,j)=Power_out;
    % j=j+1;
    % Power_l = Power_out;
    % increment = P_sumrate - P_sumrate_ini;
    % P_sumrate_ini = P_sumrate
    % end
    % P_sumrate
    % P_ga(:,:,i) = Power_out;


    gBestPosition_ga=ALPHA;
    pisition_container_ga = zeros(3,M,T,100);
    score_container_ga = zeros(100,2,100);
    for j =1:5
        [gBestScore_ga,container_ga,iter_position_ga]=ga_uav_trajectory_optimization(M,T,BETA,round(Sche),Power_out,N_0,pathLossData,gBestPosition_ga,20/sqrt(3));
        pisition_container_ga(:,:,:,5*(i-1)+j)= iter_position_ga;
        score_container_ga(:,:,5*(i-1)+j) = container_ga;
        fprintf("这是GA第%d大循环 第%d小循环",i,j);

    end

    gBestPosition_ga = iter_position_ga;
end



%%%SCA-PSO-CM Solve
gBestPosition = ALPHA;
increment_total =100;
Power_out = 0.05*ones(T,N);
Sum_ini  = 0;
i = 0;
H_out = squeeze(ALPHA(3,:,:));
%gBestPosition = pisition_container_scapsocm(:,:,:,10);
position_container_halfpower= zeros(3,M,T,100);
score_container_halfpower = zeros(100,2,100);
while i<5
    i = i+1 ;
    %Solve Scheduling Problem
    increment = 100;
    Optimal_val_ini = -100;
    A_l = 0.1*ones(M,N,T);
    while increment>0.01
    [Sche, Optimal_val,R_out,psi_A_out]=communication_scheduling_radiomap(M,N,T,A_l,Power_out,gBestPosition,BETA,pathLossData,N_0,0);
    A_l = Sche;
    increment = Optimal_val - Optimal_val_ini;
    Optimal_val_ini = Optimal_val
    end
    Optimal_val
    A_halfpower(:,:,:,i) = Sche;

    % Power_l = 0.1*ones(T,N);
    % power_container = zeros(T,N,20);
    % P_sumrate_ini = 0;
    % increment = 100;
    % j =1;
    % while increment>0.1
    % [Power_out, P_rate, P_sumrate]=Power_allocation_radiomap(M,N,T,round(Sche),gBestPosition,BETA,pathLossData,Power_l,N_0);
    % power_container(:,:,j)=Power_out;
    % j=j+1;
    % Power_l = Power_out;
    % increment = P_sumrate - P_sumrate_ini;
    % P_sumrate_ini = P_sumrate
    % end
    % P_sumrate

    %迭代求解trajectory design
    for ss = 1:5
        gBestPosition_sca = ALPHA(1:2,:,:);
        Q_sumrate_ini = -100;
        increment = 100;
        while increment>0.1
        [iter_position_sca,Q_rate,Q_sumrate] = trajectory_design(M,N,T,round(Sche),Power_out,gBestPosition_sca,H_out,BETA,L_0,N_0,20,20,d_min);
        gBestPosition_sca = iter_position_sca;
        increment = Q_sumrate - Q_sumrate_ini;
        Q_sumrate_ini = Q_sumrate
        end
        Q_sumrate
    
        %Sovle Height Problem
        H_l = squeeze(ALPHA(3,:,:));
        increment = 100;
        verticl_ini =-100;
        while increment>0.1
        [H_out,verticle_Q_rate,verticle_Q_sumrate] =vertical_trajectory_design(M,N,T,round(Sche),Power_out,gBestPosition_sca,H_l,BETA,L_0,N_0,20,20,d_min);
        H_l = H_out;
        increment = verticle_Q_sumrate - verticl_ini;
        verticl_ini = verticle_Q_sumrate
        end

    gBestPosition(1:2,:,:) = gBestPosition_sca;
    gBestPosition(3,:,:) = H_out;
    end
    for j =1:5
        [gBestScore,container,iter_position]=pso_uav_trajectory_optimization(M,T,BETA,round(Sche),Power_out,N_0,pathLossData,gBestPosition,20/sqrt(3),1);
        position_container_halfpower(:,:,:,5*(i-1)+j)= iter_position;
        score_container_halfpower(:,:,5*(i-1)+j) = container;
        gBestPosition = iter_position;
        gBestPosition = smooth_trajectory(gBestPosition, 3);
        fprintf("这是FP第%d大循环 第%d小循环",i,j);

    end

end


A_pp = zeros(M,N,T);
A_pp(1,1,1:2:T-1) = 1;
A_pp(1,2,2:2:T) = 1;
A_pp(2,3,1:2:T-1) = 1;
A_pp(2,4,2:2:T) = 1;
%%%SCA-PSO-CM Solve
gBestPosition = ALPHA;
increment_total =100;
Power_out = 0.5*ones(T,N);
Sum_ini  = 0;
i = 0;
H_out = squeeze(ALPHA(3,:,:));
%gBestPosition = pisition_container_scapsocm(:,:,:,10);
position_container_closeA= zeros(3,M,T,100);
score_container_closeA = zeros(100,2,100);
while i<5
    i = i+1 ;
    % %Solve Scheduling Problem
    % %Solve Scheduling Problem
    % increment = 100;
    % Optimal_val_ini = 0;
    % A_l = 0.1*ones(M,N,T);
    % while increment>0.01
    % [Sche, Optimal_val,R_out,psi_A_out]=communication_scheduling_radiomap(M,N,T,A_l,Power_out,gBestPosition,BETA,pathLossData,N_0,1);
    % A_l = Sche;
    % increment = Optimal_val - Optimal_val_ini;
    % Optimal_val_ini = Optimal_val
    % end
    % Optimal_val
    % A_close(:,:,:,i) = Sche;

    % Power_l = 0.1*ones(T,N);
    % power_container = zeros(T,N,20);
    % P_sumrate_ini = 0;
    % increment = 100;
    % j =1;
    % while increment>0.1
    % [Power_out, P_rate, P_sumrate]=Power_allocation_radiomap(M,N,T,round(Sche),gBestPosition,BETA,pathLossData,Power_l,N_0);
    % power_container(:,:,j)=Power_out;
    % j=j+1;
    % Power_l = Power_out;
    % increment = P_sumrate - P_sumrate_ini;
    % P_sumrate_ini = P_sumrate
    % end
    % P_sumrate
    % 
    % P_close(:,:,i) = Power_out;

    for xx =1:5
        %迭代求解trajectory design
        gBestPosition_sca = ALPHA(1:2,:,:);
        Q_sumrate_ini = -100;
        increment = 100;
        while increment>0.1
        [iter_position_sca,Q_rate,Q_sumrate] = trajectory_design(M,N,T,A_pp,Power_out,gBestPosition_sca,H_out,BETA,L_0,N_0,20,20,d_min);
        gBestPosition_sca = iter_position_sca;
        increment = Q_sumrate - Q_sumrate_ini;
        Q_sumrate_ini = Q_sumrate
        end
        Q_sumrate
    
        %Sovle Height Problem
        H_l = squeeze(ALPHA(3,:,:));
        increment = -100;
        verticl_ini = -100;
        while increment>0.1
        [H_out,verticle_Q_rate,verticle_Q_sumrate] =vertical_trajectory_design(M,N,T,A_pp,Power_out,gBestPosition_sca,H_l,BETA,L_0,N_0,20,20,d_min);
        H_l = H_out;
        increment = verticle_Q_sumrate - verticl_ini;
        verticl_ini = verticle_Q_sumrate
        end
    end
    gBestPosition(1:2,:,:) = gBestPosition_sca;
    gBestPosition(3,:,:) = H_out;

    for j =1:5
        [gBestScore,container,iter_position]=pso_uav_trajectory_optimization(M,T,BETA,A_pp,Power_out,N_0,pathLossData,gBestPosition,20/sqrt(3),1);
        position_container_closeA(:,:,:,5*(i-1)+j)= iter_position;
        score_container_closeA(:,:,5*(i-1)+j) = container;
        gBestPosition = iter_position;
        gBestPosition = smooth_trajectory(gBestPosition, 3);
        fprintf("这是RR第%d大循环 第%d小循环",i,j);

    end

end