% This script evaluates the IMU alignment algorithm on real data.

%%%%%%%%%%%%%%%%%%
%%  Load data1. %%
%%%%%%%%%%%%%%%%%%

load('data1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run IMU alignment algorithm. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[z_h1,du_tot1,P1,Psi_sb2s1] = initialization_for_smph_veh_alignment(pos1,speed1,bearing1,u1,ahrs_utc1);
[z_h1,du_tot1,Psi_sb2s1,Psi_tb2t1,index1,C_b2s1,P_tot1] = GNSSaided_smph_veh_alignment(z_h1,du_tot1,P1,Psi_sb2s1,pos1,speed1,bearing1,u1,ahrs_utc1,index_vector1,init_lla1);

[z_h2,du_tot2,P2,Psi_sb2s2] = initialization_for_smph_veh_alignment(pos2,speed2, bearing2,u2,ahrs_utc2);
[z_h2,du_tot2,Psi_sb2s2,Psi_tb2t2,index2,C_b2s2,P_tot2] = GNSSaided_smph_veh_alignment(z_h2,du_tot2,P2,Psi_sb2s2,pos2,speed2,bearing2,u2,ahrs_utc2,index_vector2,init_lla1);

[z_h3,du_tot3,P3,Psi_sb2s3] = initialization_for_smph_veh_alignment(pos3,speed3,bearing3,u3,ahrs_utc3);
[z_h3,du_tot3,Psi_sb2s3,Psi_tb2t3,index3,C_b2s3,P_tot3]= GNSSaided_smph_veh_alignment(z_h3,du_tot3,P3,Psi_sb2s3,pos3,speed3,bearing3,u3,ahrs_utc3,index_vector3,init_lla3);

%%%%%%%%%%%%%%%%%%%%%%&%%%%%%%%%%%%%%%%
%%  Run standard INS for smartphone. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[z_h1s,du_tot1s,P1s] = initialization_for_standard_INS(pos1,speed1,bearing1,u1,ahrs_utc1);
[z_h1s,du_tot1s,P1s,index1s] = standard_GNSSaided_INS(z_h1s,du_tot1s,P1s,pos1,speed1,bearing1,u1,ahrs_utc1,index_vector1,init_lla1);

[z_h2s,du_tot2s,P2s] = initialization_for_standard_INS(pos2,speed2,bearing2,u2,ahrs_utc2);
[z_h2s,du_tot2s,P2s,index2s] = standard_GNSSaided_INS(z_h2s,du_tot2s,P2s,pos2,speed2,bearing2,u2,ahrs_utc2,index_vector2,init_lla2);

[z_h3s,du_tot3s,P3s] = initialization_for_standard_INS(pos3,speed3,bearing3,u3,ahrs_utc3);
[z_h3s,du_tot3s,P3s,index3s] = standard_GNSSaided_INS(z_h3s,du_tot3s,P3s,pos3,speed3,bearing3,u3,ahrs_utc3,index_vector3,init_lla3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate errors in Euler angle estimates. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectors of standard errors.
Euler_se_1 = zeros(3,length(Psi_tb2t1)); 
Euler_se_2 = zeros(3,length(Psi_tb2t2));
Euler_se_3 = zeros(3,length(Psi_tb2t3));

for n = 1:length(Psi_tb2t1)
    [~,ind] = min(abs(ahrs_utc-ahrs_utc1(n)));
    temp_euler = [ref_z_h(7:9,ind)-Psi_tb2t1(:,n) 2*pi*ones(3,1)+ref_z_h(7:9,ind)-Psi_tb2t1(:,n) -2*pi*ones(3,1)+ref_z_h(7:9,ind)-Psi_tb2t1(:,n)];
    [~,ind2] = min(abs(temp_euler'));
    Euler_se_1(:,n) = [temp_euler(1,ind2(1)) temp_euler(2,ind2(2)) temp_euler(3,ind2(3))]';
    
    [~,ind] = min(abs(ahrs_utc-ahrs_utc2(n)));
    temp_euler = [ref_z_h(7:9,ind)-Psi_tb2t2(:,n) 2*pi*ones(3,1)+ref_z_h(7:9,ind)-Psi_tb2t2(:,n) -2*pi*ones(3,1)+ref_z_h(7:9,ind)-Psi_tb2t2(:,n)];
    [~,ind2] = min(abs(temp_euler'));
    Euler_se_2(:,n) = [temp_euler(1,ind2(1)) temp_euler(2,ind2(2)) temp_euler(3,ind2(3))]';
    
    [~,ind] = min(abs(ahrs_utc-ahrs_utc3(n)));
    temp_euler = [ref_z_h(7:9,ind)-Psi_tb2t3(:,n) 2*pi*ones(3,1)+ref_z_h(7:9,ind)-Psi_tb2t3(:,n) -2*pi*ones(3,1)+ref_z_h(7:9,ind)-Psi_tb2t3(:,n)];
    [~,ind2] = min(abs(temp_euler'));
    Euler_se_3(:,n) = [temp_euler(1,ind2(1)) temp_euler(2,ind2(2)) temp_euler(3,ind2(3))]';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Euler angle errors. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sqrt(mean(Euler_se_1.^2,2))*180/pi) % [degree]
disp(sqrt(mean(Euler_se_2.^2,2))*180/pi) % [degree]
disp(sqrt(mean(Euler_se_3.^2,2))*180/pi) % [degree]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate position drifts during GNSS outage. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pos_error1,pos_error_standard1,error_time1] = pos_drift_at_gnss_outtage(pos1,speed1,bearing1,ahrs_utc1,u1,index_vector1,init_lla1,ref_z_h,ahrs_utc,index1,index1s);
[pos_error2,pos_error_standard2,error_time2] = pos_drift_at_gnss_outtage(pos2,speed2,bearing2,ahrs_utc2,u2,index_vector2,init_lla2,ref_z_h,ahrs_utc,index2,index2s);
[pos_error3,pos_error_standard3,error_time3] = pos_drift_at_gnss_outtage(pos3,speed3,bearing3,ahrs_utc3,u3,index_vector3,init_lla3,ref_z_h,ahrs_utc,index3,index3s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plot position drifts during GNSS outage. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_drift_at_gnss_outtage(pos_error1,pos_error_standard1,error_time1,pos_error2,pos_error_standard2,error_time2,pos_error3,pos_error_standard3,error_time3)

%%%%%%%%%%%%%%%%%%%%%
%% Plot trajectory %%
%%%%%%%%%%%%%%%%%%%%%

close all
figure('units','normalized','position',[.1 .1 0.3 .22]);
clf
plot(ref_z_h(2,:),ref_z_h(1,:),'.')
hold on
title('title')
xlabel('xlabel')
ylabel('ylabel')
axis equal
% axis([0 1200 0 500])

%%

close all
fig1 = figure('units','normalized','position',[.1 .1 0.3 .15]);
plot(ahrs_utc-ahrs_utc(1),sqrt(sum(ref_z_h(4:5,:).^2,1))*3.6,'LineWidth',2)
title('title')
xlabel('xlabel')
ylabel('ylabel')

