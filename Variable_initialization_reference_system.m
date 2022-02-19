% Initialization file for the GPS aided inertial navigation system when 
% processing real world data.
%
% Note: All setting are stored in the global variable “simdata” to be 
% accessible from all functions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              GENERAL PARAMETERS         %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global simdata;

simdata.earth_rotation = 7.292115*1e-5;                                         % See pages 44,45 in Groves (2008).
simdata.Omega_ie2e = Skew_symmetric_matrix([0 0 simdata.earth_rotation]);       % See eq. (5.18) and (A.27) in Groves(2008).                                                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             FILTER PARAMETERS           %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Process noise covariance (Q)
simdata.sigma_acc = 10*9.8*80*10^-6*ones(3,1);          % [m/(s^2*sqrt(Hz))], models the measurement noise of the accelerometer. See product sheet for 3DM-GX3 35 (a factor 10 has been added).
simdata.sigma_gyro = 10*0.03*ones(3,1)*pi/180;          % [rad/(s*sqrt(Hz))], models the measurement noise of the gyroscope. See product sheet for 3DM-GX3 35 (a factor 10 has been added).
simdata.sigma_acc_bias=1e-8;                            % [m/(s^3*sqrt(Hz))], models how the bias in the accelerometer changes. 
simdata.sigma_gyro_bias=1e-8*pi/180;                    % [rad/(s^2*sqrt(Hz))], models how the bias in the gyroscope changes. 

% Measurement noise covariance (R) 
% Position error std
simdata.sigma_gps_pos_horiz=3/sqrt(2*log(20));          % [m]      
simdata.sigma_gps_pos_vert=3/sqrt(2*log(20));           % [m]      
simdata.sigma_gps_speed=0.3;                            % [m/s]     
simdata.sigma_vel_pseudo_obs_side=1000000;              % [m/s]
simdata.sigma_vel_pseudo_obs_down=1000000;              % [m/s]

% Initial uncertainties (standard deviations)  
simdata.init_sd(1)=3/sqrt(3);                           % Position [m]
simdata.init_sd(2)=3/sqrt(3);                           % Velocity [m/s]
simdata.init_sd(3:5)=(pi/180*[4 4 10]');                % Attitude (roll,pitch,yaw) [rad]
simdata.init_sd(6)=0.05;                                % Accelerometer biases [m/s^2]
simdata.init_sd(7)=(0.5*pi/180);                        % Gyro biases [rad/s]                               


