% This function calculates the position drifts at GNSS outages.
%
% Inputs:  
%               - pos1: Position measurements from smartphone.
%               - speed1: Speed measurements from smartphone.
%               - bearing1: Bearing measurements from smartphone.
%               - ahrs_utc1: Sampling instances from smartphone.
%               - u1: IMU measurements from smartphone.
%               - index_vector1: Specifying the mapping between IMU and GNSS data.
%               - init_lla: Lat-Long-Alt in origin of tangent frame.
%               - ref_z_h: Navigation state from reference system.
%               - ahrs_utc: Sampling instances from reference system. 
%               - index1: Identified initialization particle for IMU alignment.
%               - index1s: Identified initialization particle for standard GNSS-aided INS.
%
% Outputs: 
%               - pos_error: Position error from IMU alignment method.
%               - pos_error_standard: Position error from standard GNSS-aided INS.
%               - error_time: Time since start of outage.
%

function [pos_error,pos_error_standard,error_time] = pos_drift_at_gnss_outtage(pos1,speed1,bearing1,ahrs_utc1,u1,index_vector1,init_lla1,ref_z_h,ahrs_utc,index1,index1s)

% Time at which the outages are to start.
start_time = 4000;

% Length of each outage in sampling instances.
stop_length = 1020;
iterations = floor((size(u1,2)-(start_time+stop_length))/stop_length);

error_time = zeros(iterations,stop_length); 
pos_error = zeros(iterations,stop_length);
pos_error_standard = zeros(iterations,stop_length);

cumsum_index_vector = zeros(size(index_vector1));
cumsum_index_vector(1) = index_vector1(1);
for n = 2:length(index_vector1)
    cumsum_index_vector(n) = cumsum_index_vector(n-1)+index_vector1(n);
end

% Iterate over each simulated outage.
for iter = 1:iterations
 
disp(iter)

% Time for when outage starts and stops.
ahrs_t_1 = start_time+stop_length*iter;
ahrs_t_2 = ahrs_t_1 + stop_length;

stop_index = sum(index_vector1(1:ahrs_t_1));

ref_pos_out = pos1(:,1:stop_index);
ref_speed_out = speed1(1:stop_index);
ref_bearing_out = bearing1(1:stop_index);

u_out = u1(:,1:ahrs_t_2);
index_vector_out = index_vector1;
index_vector_out(ahrs_t_1+1:end) = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run IMU alignment method. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[z_h1,du_tot1,P1,Psi_sb2s1] = initialization_for_smph_veh_alignment(ref_pos_out,ref_speed_out,ref_bearing_out,u_out,ahrs_utc1);
z_h1 = z_h1(:,:,index1);
du_tot1 = du_tot1(:,:,index1);
P1 = P1(:,:,index1);
Psi_sb2s1 = Psi_sb2s1(:,:,index1);
[ref_z_h_out2,~,~,~,~,~] = GNSSaided_smph_veh_alignment(z_h1,du_tot1,P1,Psi_sb2s1,ref_pos_out,ref_speed_out,ref_bearing_out,u_out,ahrs_utc1,index_vector_out,init_lla1);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run standard GNSS-aided INS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[z_h1s,du_tot1s,P1s] = initialization_for_standard_INS(ref_pos_out,ref_speed_out,ref_bearing_out,u_out,ahrs_utc1);
z_h1s = z_h1s(:,:,index1s);
du_tot1s = du_tot1s(:,:,index1s);
P1s = P1s(:,:,index1s);  
[ref_z_h_out12,~,~,~] = standard_GNSSaided_INS(z_h1s,du_tot1s,P1s,ref_pos_out,ref_speed_out,ref_bearing_out,u_out,ahrs_utc1,index_vector_out,init_lla1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save position drifts during outage %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = ahrs_t_1:ahrs_t_2
    [~,ind] = min(abs(ahrs_utc-ahrs_utc1(j)));
    error_time(iter,j-ahrs_t_1+1) = ahrs_utc1(j)-ahrs_utc1(ahrs_t_1);
    pos_error(iter,j-ahrs_t_1+1) = sqrt(sum((ref_z_h(1:2,ind)-ref_z_h_out2(1:2,j)).^2));
    pos_error_standard(iter,j-ahrs_t_1+1) = sqrt(sum((ref_z_h(1:2,ind)-ref_z_h_out12(1:2,j)).^2));
end

end

end

