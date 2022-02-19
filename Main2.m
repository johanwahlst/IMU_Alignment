% This script evaluates the IMU alignment algorithm on simulated data.

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load simulated data %%
%%%%%%%%%%%%%%%%%%%%%%%%%

load('simdata')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate noise and run IMU Alignment. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Variable_initialization_smartphone

Numb_of_simulations = 5;        % Number of simulations.
dT_GNSS_fact = 20;              % Factor of IMU and GNSS sample rate. 

% Angle error variables.
Psi_sb2s_errors = zeros(3,size(ref_z_h,2),Numb_of_simulations);
Psi_tb2t_errors = zeros(3,size(ref_z_h,2),Numb_of_simulations);
Psi_ts2t_errors = zeros(3,size(ref_z_h,2),Numb_of_simulations);

Cov_tb2t_tot = zeros(3,3,size(ref_z_h,2),Numb_of_simulations);
Cov_sb2s_tot = zeros(3,3,size(ref_z_h,2),Numb_of_simulations);

for n = 1:Numb_of_simulations
    
    disp('iteration')
    disp(n)
    
    Psi_sb2s = [2*pi*(rand-1/2) pi*(rand-1/2) 2*pi*(rand-1/2)]';
    C_s2b = Rot_Mat_Fnc(Psi_sb2s);
    
    Psi_ts2t = zeros(3,size(ref_z_h,2));

    for m = 1:size(ref_z_h,2)
        C_t2b = Rot_Mat_Fnc(ref_z_h(7:9,m));
        C_s2t = C_t2b'*C_s2b; 
        Psi_ts2t(1,m)=atan2(C_s2t(3,2),C_s2t(3,3));                       % roll.
        Psi_ts2t(2,m)=-asin(C_s2t(3,1));                                  % pitch. 
        Psi_ts2t(3,m)=atan2(C_s2t(2,1),C_s2t(1,1));  
    end
    
    % Add noise GNSS measurements.. 
    pos1 = ref_z_h(1:3,1:dT_GNSS_fact:end)+normrnd(0,simdata.sigma_gps_pos_horiz,[3 length(1:dT_GNSS_fact:length(vel))]);
    pos1(3,:) = pos1(3,:)+normrnd(0,0.2,[1 size(pos1,2)]);
    speed1 = sqrt(sum(ref_z_h(4:5,1:dT_GNSS_fact:end).^2,1))+normrnd(0,simdata.sigma_gps_speed,[1 size(pos1,2)]);
    bearing1 = atan2(ref_z_h(5,1:dT_GNSS_fact:end),ref_z_h(4,1:dT_GNSS_fact:end));
    for m = 1:size(bearing1,2)
        bearing1(m) = bearing1(m)+normrnd(0,simdata.sigma_gps_speed/speed1(m));
    end
    
    % Add noise to IMU measurements.
    u_bias = [(rand(3,1)-1/2); 0.003*(2*rand(3,1)-1)];
    u1 = u+repmat(u_bias,[1 size(u,2)]);
    u_noise = zeros(size(u1));
    u_noise(1:3,:) = normrnd(0,simdata.sigma_acc(1),[3 size(u,2)]);
    u_noise(4:6,:) = normrnd(0,simdata.sigma_gyro(1),[3 size(u,2)]);
    u1 = u1 + u_noise; 
    u1(1:3,:) = C_s2b'*u1(1:3,:);
    u1(4:6,:) = C_s2b'*u1(4:6,:);
    
    index_vector1 = zeros(1,length(ahrs_utc1));
    index_vector1(1:dT_GNSS_fact:end) = 1; 
 
    % Run IMU Alignment.
    [z_h1,du_tot1,P1,Psi_sb2s1] = initialization_for_smph_veh_alignment(pos1,speed1,bearing1,u1,ahrs_utc1);
    [z_h1,du_tot1,Psi_sb2s1,Psi_tb2t1,index1,C_b2s1,P_tot1] = GNSSaided_smph_veh_alignment(z_h1,du_tot1,P1,Psi_sb2s1,pos1,speed1,bearing1,u1,ahrs_utc1,index_vector1,init_lla);
    
    % Save angle errors.
    Psi_tb2t_errors(:,:,n) = Psi_tb2t1-Psi_tb2t;
    Psi_ts2t_errors(:,:,n) = z_h1(7:9,:,index1)-Psi_ts2t;
    for m = 1:size(pos,2)
        Error_rot_mat = C_b2s1(:,:,m)*C_s2b;
        Psi_sb2s_errors(1,m,n)=atan2(Error_rot_mat(3,2),Error_rot_mat(3,3));                       
        Psi_sb2s_errors(2,m,n)=-asin(Error_rot_mat(3,1));                                   
        Psi_sb2s_errors(3,m,n)=atan2(Error_rot_mat(2,1),Error_rot_mat(1,1));                       
    end
    
    % Save covariance of angle error.
    Cov_tb2t = zeros(3,3,size(ref_z_h,2));
    Cov_sb2s = zeros(3,3,size(ref_z_h,2));
    for m = 1:size(ref_z_h,2)
        Cov_ts2t = P_tot1(7:9,7:9,m);
        Cov_sb2s_ts2t = P_tot1(16:18,7:9,m);
        Cov_sb2s(:,:,m) = P_tot1(16:18,16:18,m);
        C_s2t = Rot_Mat_Fnc(z_h1(7:9,m,index1))';
        Cov_tb2t(:,:,m) = Cov_ts2t+C_s2t*Cov_sb2s_ts2t+(C_s2t*Cov_sb2s_ts2t)'+C_s2t*Cov_sb2s(:,:,m)*C_s2t'; 
    end
    
    Cov_tb2t_tot(:,:,:,n) = Cov_tb2t;    
    Cov_sb2s_tot(:,:,:,n) = Cov_sb2s;    
    
end

% Only allow -pi < angle errors < pi.
Psi_tb2t_errors(Psi_tb2t_errors>pi) = Psi_tb2t_errors(Psi_tb2t_errors>pi)-2*pi;
Psi_tb2t_errors(Psi_tb2t_errors<-pi) = Psi_tb2t_errors(Psi_tb2t_errors<-pi)+2*pi;

Psi_sb2s_errors(Psi_sb2s_errors>pi) = Psi_sb2s_errors(Psi_sb2s_errors>pi)-2*pi;
Psi_sb2s_errors(Psi_sb2s_errors<-pi) = Psi_sb2s_errors(Psi_sb2s_errors<-pi)+2*pi;

Psi_ts2t_errors(Psi_ts2t_errors>pi) = Psi_ts2t_errors(Psi_ts2t_errors>pi)-2*pi;
Psi_ts2t_errors(Psi_ts2t_errors<-pi) = Psi_ts2t_errors(Psi_ts2t_errors<-pi)+2*pi;

%%%%%%%%%%%%
%% Plot 1 %%
%%%%%%%%%%%%

close all
figure('units','normalized','position',[.1 .1 0.3 .22]);
clf
plot(pos(1,:),pos(2,:)+445,'.')
hold on
title('title')
xlabel('xlabel')
ylabel('ylabel')
axis equal
axis([0 1200 0 500])
annotation('textbox',[0.48 0.77 0.15 0.06],'String','textlabel1','Edgecolor',[1 1 1])
annotation('textbox',[0.7 0.2 0.15 0.06],'String','textlabel2','Edgecolor',[1 1 1])
annotation('textbox',[0.74 0.5 0.15 0.06],'String','textlabel3','Edgecolor',[1 1 1])
annotation('textbox',[0.2 0.4 0.2 0.095],'String',{'textlabel4','textlabel5'},'Edgecolor',[1 1 1])

%%%%%%%%%%%%
%% Plot 2 %%
%%%%%%%%%%%%

close all
fig1 = figure('units','normalized','position',[.1 .1 0.3 .15]);
plot(0:dt:dt*(length(vel)-1),sqrt(sum(vel(1:2,:).^2,1))*3.6,'LineWidth',2)
title('title')
xlabel('xlabel')
ylabel('ylabel')

%%%%%%%%%%%%
%% Plot 3 %%
%%%%%%%%%%%%

plot_angle_errors_tb2t(Psi_tb2t_errors,Cov_tb2t_tot)
