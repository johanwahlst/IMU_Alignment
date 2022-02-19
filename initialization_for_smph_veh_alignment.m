% This function initializes variables used in the IMU alignment algorithm. 
%
% Inputs: 
%               - pos: position measurements.
%               - speed: speed measurements.
%               - bearing: bearing measurements.
%               - u: IMU measurements (acc + gyro).
%               - ahrs_utc: sampling instances.
%
% Outputs: 
%               - z_h: navigation state.
%               - du_tot: Estimated IMU bias.
%               - P: state covariance.
%               - Psi_sb2s: Euler angles describign the smph-to-vehicle orientation.

function [z_h,du_tot,P,Psi_sb2s,grid_points] = initialization_for_smph_veh_alignment(pos,speed,bearing,u,ahrs_utc)

Variable_initialization_smartphone

            % Number of grid points in MPF for: 
n1 = 1;     % Smartphone roll.
n2 = 1;     % Smartphone pitch.
n3 = 8;     % Smartphone yaw.

Numb_of_particles = n1*n2*n3; 
grid_points = [n1 n2 n3];

z_h = zeros(9,length(u),Numb_of_particles);

mean_force=mean(u(1:3,ahrs_utc<ahrs_utc(1)+1),2);
init_roll_pitch = [atan2(-mean_force(2),-mean_force(3)) atan2(mean_force(1),norm(mean_force(2:3)))];

% Initialize navigation state.
z_h(1:3,1,:) = repmat(pos(:,1),[1 Numb_of_particles]);                                       % Initialize position.
z_h(4:5,1,:) = repmat(speed(1)*[cos(bearing(1)) sin(bearing(1))]',[1 Numb_of_particles]);  % Initialize velocity.
z_h(7:9,1,:) = attitude_init(n1,n2,n3,init_roll_pitch(1),init_roll_pitch(2));                 % Initialize smartphone attitude.

du_tot = zeros(6,length(u),Numb_of_particles);

% Initialize state covariance matrix.
P=zeros(18,18,Numb_of_particles);

P(1:3,1:3,:)=repmat(simdata.init_sd(1)^2*eye(3),[1 1 Numb_of_particles]);
P(4:6,4:6,:)=repmat(simdata.init_sd(2)^2*eye(3),[1 1 Numb_of_particles]);
P(7:9,7:9,:)=repmat(diag(simdata.init_sd(3:5)).^2,[1 1 Numb_of_particles]);
P(10:12,10:12,:)=repmat(simdata.init_sd(6)^2*eye(3),[1 1 Numb_of_particles]);
P(13:15,13:15,:)=repmat(simdata.init_sd(7)^2*eye(3),[1 1 Numb_of_particles]);
% P(16:18,16:18,:)=repmat(simdata.init_sd(8)^2*eye(3),[1 1 Numb_of_particles]);
for particle = 1:Numb_of_particles
    P(16:18,16:18,particle)=Rot_Mat_Fnc(z_h(7:9,1,particle))*diag(simdata.init_sd(3:5)+simdata.init_sd(8:10)).^2*Rot_Mat_Fnc(z_h(7:9,1,particle))';
end

% Initialize smartphone-to-vehicle orientation.
Psi_sb2s = zeros(3,min((length(u)-1),simdata.Numb_of_samples_at_init),Numb_of_particles); % Vector for smartphone-vehicle attitude.

for particle = 1:Numb_of_particles
    C_t2s = Rot_Mat_Fnc(z_h(7:9,1,particle));
    C_t2b = Rot_Mat_Fnc([0 0 bearing(1)]);
    C_b2s = C_t2s*C_t2b';
    Psi_sb2s(1,1,particle)=atan2(C_b2s(3,2),C_b2s(3,3));
    Psi_sb2s(2,1,particle)=-asin(C_b2s(3,1));
    Psi_sb2s(3,1,particle)=atan2(C_b2s(2,1),C_b2s(1,1));
end

end

