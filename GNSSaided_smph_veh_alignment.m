% This is the Main function for the IMU alignment algorithm.
%
% Inputs:
%           - z_h: smartphone state (pos, vel, att).
%           - du_tot: IMU bias (acc + gyro).
%           - P: State covariance matrix.
%           - Psi_sb2s: Euler angles describing the smph-to-veh orientation.
%           - pos: Postiion measurements.
%           - speed: Speed measurements.
%           - bearing: Bearing measurements.
%           - u: IMU measurements (acc + gyro).
%           - ahrs_utc: Sampling instances.
%           - index_vector: Specifying the mapping between IMU and GNSS data.
%           - init_lla: Lat-Long-Alt in origin of tangent frame.
%
% Outputs: 
%           - z_h: Navigation solution.
%           - du_tot: Estimated IMU bias.
%           - Psi_sb2s: Euler angles describing the smph-to-veh orientation.
%           - Psi_tb2t: Euler angles describing the tangent-to-veh orientation.
%           - index: Identified initialization particle.
%           - du_out: Estimated IMU bias.
%           - C_b2s: Rotation matrix.
%           - P_tot: State covariance matrix.

function [z_h,du_tot,Psi_sb2s,Psi_tb2t,index,C_b2s,P_tot] = GNSSaided_smph_veh_alignment(z_h,du_tot,P,Psi_sb2s,pos,speed,bearing,u,ahrs_utc,index_vector,init_lla)

%%%%%%%%%%%%%%%%%%%%
%% Initialization %%
%%%%%%%%%%%%%%%%%%%%

global simdata

Numb_of_particles = size(z_h,3);

T_tot = zeros(1,Numb_of_particles);

eps = zeros(8,Numb_of_particles); % Normalized innovations.

C_b2s = zeros(3,3,(min((length(u)-1),simdata.Numb_of_samples_at_init)+1),Numb_of_particles);
for n = 1:Numb_of_particles
    C_b2s(:,:,1,n) = Rot_Mat_Fnc(Psi_sb2s(:,1,n))';
end

m = 2;

P_tot = zeros(18,18,(min((length(u)-1),simdata.Numb_of_samples_at_init)+1),Numb_of_particles);
for particle = 1:Numb_of_particles
    P_tot(:,:,1,particle) = P(:,:,particle);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run INS for smartphone vehicle alignment in MPF %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reverseStr = '';
for k = 1:min((length(u)-1),simdata.Numb_of_samples_at_init)
    y = zeros(9,1);
    if(index_vector(k+1))
        y = [pos(:,m); speed(m); bearing(m)]; 
        m = m + 1;
    end
    if mod(k,20)==0
        msg = sprintf('Processed %d/%d', k, length(u)-1);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    for particle = 1:Numb_of_particles
        [z_h(:,k+1,particle),~,P(:,:,particle),~,du_tot(:,k+1,particle),eps(:,particle),C_b2s(:,:,k+1,particle)] = GNSSaidedINS_one_iteration_18_state(y,z_h(:,k,particle),u(:,k),P(:,:,particle),index_vector(k+1),ahrs_utc(k+1)-ahrs_utc(k),init_lla,du_tot(:,k,particle),C_b2s(:,:,k,particle));
        P_tot(:,:,k+1,particle) = P(:,:,particle);
    end
    sigma_T = 1;
    T_tot_temp = zeros(1,Numb_of_particles);
    while(sum(T_tot_temp(:)==0)) 
        for particle = 1:Numb_of_particles
            T_tot_temp(particle) = T_tot_temp(particle) + log(prod(normpdf(eps(:,particle),0,sigma_T)));
        end
        sigma_T = sigma_T*4;
    end
    T_tot = T_tot+T_tot_temp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the initial attitude. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,index] = max(T_tot);
P_tot_temp = zeros(18,18,length(u));
P_tot_temp(:,:,1:min((length(u)-1),simdata.Numb_of_samples_at_init)+1) = P_tot(:,:,:,index);
P_tot = P_tot_temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Continue INS for smartphone with the identified initial attitude. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reduce the vector of Euler angles.
Psi_sb2s = Psi_sb2s(:,:,index);
Psi_sb2s = [Psi_sb2s zeros(3,(length(u)-1)-(min((length(u)-1),simdata.Numb_of_samples_at_init)+1))];

C_b2s_temp = zeros(3,3,size(z_h,2));
C_b2s_temp(:,:,1:(min((length(u)-1),simdata.Numb_of_samples_at_init)+1)) = C_b2s(:,:,:,index);
C_b2s = C_b2s_temp;

for k = (min((length(u)-1),simdata.Numb_of_samples_at_init)+1):(length(u)-1)
    y = zeros(9,1);
    if(index_vector(k+1))
        y = [pos(:,m); speed(m); bearing(m)]; 
        m = m + 1;
    end
    if mod(k,20)==0
        msg = sprintf('Processed %d/%d', k, length(u)-1);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    for particle = index
        [z_h(:,k+1,particle),~,P(:,:,particle),~,du_tot(:,k+1,particle),eps(:,particle),C_b2s(:,:,k+1)] = GNSSaidedINS_one_iteration_18_state(y,z_h(:,k,particle),u(:,k),P(:,:,particle),index_vector(k+1),ahrs_utc(k+1)-ahrs_utc(k),init_lla,du_tot(:,k,particle),C_b2s(:,:,k));
        P_tot(:,:,k+1) = P(:,:,particle);
    end
end

% Calculate Euler angles.
Psi_tb2t = zeros(3,length(u));
for k = 1:size(u,2)
    C_t2s = Rot_Mat_Fnc(z_h(7:9,k,index));
    C_b2s_temp = C_b2s(:,:,k);
    C_s2b = C_b2s_temp';
    Psi_sb2s(1,k)=atan2(C_b2s_temp(3,2),C_b2s_temp(3,3));                      
    Psi_sb2s(2,k)=-asin(C_b2s_temp(3,1));                                  
    Psi_sb2s(3,k)=atan2(C_b2s_temp(2,1),C_b2s_temp(1,1));    
    C_b2t = C_t2s'*C_s2b';
    Psi_tb2t(1,k)=atan2(C_b2t(3,2),C_b2t(3,3));                       
    Psi_tb2t(2,k)=-asin(C_b2t(3,1));   
    Psi_tb2t(3,k)=atan2(C_b2t(2,1),C_b2t(1,1));  
end

fprintf('\n')

end

