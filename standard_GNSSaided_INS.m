% This function integrates GNSS and IMU data in a GNSS-aided INS.
% 
% Inputs:
%           - z_h: navigation state (pos, vel, att).
%           - du_in: IMU bias (acc + gyro).
%           - P: State covariance matrix.
%           - pos1: Position measurements.
%           - speed1: Speed measurements.
%           - bearing1: Bearing measurements.
%           - u1: IMU measurements.
%           - ahrs_utc1: Sampling instances.
%           - index_vector: Specifying the mapping between IMU and GNSS data.
%           - init_lla: Lat-Long-Alt in origin of tangent frame.
%
% Outputs: 
%           - z_h: Navigation solution.
%           - du_in: IMU bias.
%           - P: State covariance matrix.
%           - index: The identified initialization particle.

function [z_h,du_in,P,index] = standard_GNSSaided_INS(z_h,du_in,P,pos1,speed1,bearing1,u1,ahrs_utc1,index_vector1,init_lla1)

global simdata

Numb_of_particles = size(z_h,3); 

T_tot = zeros(1,Numb_of_particles);

ref_pred_err = zeros(5,length(u1),Numb_of_particles);
eps = zeros(5,Numb_of_particles);

m = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Running standard INS for smartphone data. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reverseStr = '';
for k = 1:min(length(u1)-1,simdata.Numb_of_samples_at_init)
    y = zeros(5,1);
    if(index_vector1(k+1))
        y = [pos1(:,m); speed1(m); bearing1(m)]; 
        eps = zeros(5,Numb_of_particles);
        m = m + 1; 
    end
    if mod(k,20)==0
        msg = sprintf('Processed %d/%d', k, length(u1)-1);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    for particle = 1:Numb_of_particles
        [z_h(:,k+1,particle),~,P(:,:,particle),ref_pred_err(:,k+1,particle),du_in(:,k+1,particle),eps(:,particle)] = GNSSaidedINS_one_iteration_15_state(y,u1(:,k),index_vector1(k+1),ahrs_utc1(k+1)-ahrs_utc1(k),init_lla1,P(:,:,particle),z_h(:,k,particle),du_in(:,k,particle));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Continue INS with the identified initial attitude. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = (min(length(u1)-1,simdata.Numb_of_samples_at_init)+1):(length(u1)-1)
    y = zeros(5,1);
    if(index_vector1(k+1))
        y = [pos1(:,m); speed1(m); bearing1(m)]; 
        eps = zeros(5,Numb_of_particles);
        m = m + 1;
    end
    if mod(k,20)==0
        msg = sprintf('Processed %d/%d', k, length(u1)-1);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    for particle = index
        [z_h(:,k+1,particle),~,P(:,:,particle),ref_pred_err(:,k+1,particle),du_in(:,k+1,particle),eps(:,particle)] = GNSSaidedINS_one_iteration_15_state(y,u1(:,k),index_vector1(k+1),ahrs_utc1(k+1)-ahrs_utc1(k),init_lla1,P(:,:,particle),z_h(:,k,particle),du_in(:,k,particle));
    end
end

fprintf('\n')

end

