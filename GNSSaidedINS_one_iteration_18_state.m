% This function performs one iteration of a GNSS-aided INS with IMU alignment.
% 
% Inputs:
%           - y: GNSS measurements (pos, vel, bear).
%           - z_in: smartphone state (pos, vel, att).
%           - u: IMU measurements (acc + gyro).
%           - P_in: State covariance matrix.
%           - index_vector: Specifying the mapping between IMU and GNSS data.
%           - dT: Time between samples. 
%           - init_lla: Lat-Long-Alt in origin of tangent frame.
%           - du: Estimated IMU bias.
%           - C_b2s: Rotation matrix.
%
% Outputs: 
%           - z_h: Navigation solution.
%           - u_h: IMU measurement with estimated bias subtracted.
%           - P_out: State covariance matrix.
%           - pred_err: Prediction error.
%           - du_out: Estimated IMU bias.
%           - eps: Normalized innovations.
%           - C_b2s: Rotation matrix.

function [z_h,u_h,P_out,pred_err,du_out,eps,C_b2s]=GNSSaidedINS_one_iteration_18_state(y,z_in,u,P_in,index_vector,dT,init_lla,du,C_b2s)

C_s2b = C_b2s';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Run filter                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = u-du; 

% Navigation equation update.
[z_in]=propagate_navigation_estimates(z_in,u,init_lla,dT);
% Calculate filter matrices.
[F,G,Q]=calculate_filter_matrices(z_in,u,init_lla,dT);

% Time update of estimated variance P(k+1|k), P=F*P*F'+G*Q*G';
P_in=F*P_in*F'+G*Q*G'*dT; % See eq. (4.111) in Farrell (2008) for the approximation giving the second term.

% A standard measurement update is performed if there are GPS measurements.
if(index_vector)
    
   [K,H,R]=Gain(P_in,z_in,C_b2s);                                                       % Calculation of gain matrix. 
   bearing_h = atan2(z_in(5),z_in(4));                                                  % Calculating estimated bearing. 

   % Adjusting the measured and estimated bearing to lie in the same interval.
   if(abs(bearing_h-(y(5)-2*pi))<abs(bearing_h-y(5)))           
       y(5) = y(5)-2*pi;
   elseif(abs(bearing_h-(y(5)+2*pi))<abs(bearing_h-y(5)))       
       y(5) = y(5)+2*pi;
   end
   
   % Calculate pseudo observations.
   C_t2s = Rot_Mat_Fnc(z_in(7:9));
   C_t2b =  C_s2b*C_t2s;
   pred_err_pseudo_obs = [zeros(2,1) eye(2)]*C_s2b*C_t2s*z_in(4:6);
   
   pred_err=[z_in(1:3); norm(z_in(4:5)); bearing_h; pred_err_pseudo_obs; atan2(C_t2b(2,3),C_t2b(3,3))]-[y(:); zeros(3,1)];   % Calculation of prediction error.
   dx=[zeros(9,1); du; zeros(3,1)]+K*pred_err;                                          % Error estimation: delta_x(k) = [delta_z(k) delta_u(k)].

   dz=dx(1:9);                                                                          % Navigation errors.
   du_out=dx(10:15);                                                                    % Sensor bias.
   dPsi=dx(16:18);
   
   eps = (sqrtm(H*P_in*H'+R))\pred_err;                                                 % Normalized innovations.
   P_out=(eye(18)-K*H)*P_in;                                                            % Update covariance matrix P(k|k).
   [z_h,C_b2s]=correct_navigation_estimates(z_in,dz,dPsi,C_b2s);                        % Error compensation.
   
else
    
   [K,H,R]=Gain2(P_in,z_in,C_b2s); 
   C_t2s = Rot_Mat_Fnc(z_in(7:9));
   C_t2b =  C_s2b*C_t2s;
   pred_err_pseudo_obs = [zeros(2,1) eye(2)]*C_s2b*C_t2s*z_in(4:6);
   pred_err_pseudo_obs = [pred_err_pseudo_obs; atan2(C_t2b(2,3),C_t2b(3,3))];
   pred_err=pred_err_pseudo_obs;   % Calculation of prediction error.
   dx=[zeros(9,1); du; zeros(3,1)]+K*pred_err;                                          % Error estimation: delta_x(k) = [delta_z(k) delta_u(k)].

   dz=dx(1:9);                                                                          % Navigation errors.
   du_out=dx(10:15);                                                                    % Sensor bias.
   dPsi=dx(16:18);
   
   eps = [zeros(5,1); (sqrtm(H*P_in*H'+R))\pred_err];                                   % Normalized innovations.
   P_out=(eye(18)-K*H)*P_in;                                                            % Update covariance matrix P(k|k).
   [z_h,C_b2s]=correct_navigation_estimates(z_in,dz,dPsi,C_b2s);                        % Error compensation.
   pred_err = [zeros(5,1); pred_err];                                                   % Calculation of prediction error.
   
end

% Compensate for sensor bias.
u_h = u-du_out;                                                                         % See page 152 in Groves (2008).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          SUB-FUNCTIONS                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate_filter_matrices  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,G,Q]=calculate_filter_matrices(z_h,u,init_lla,dT)

global simdata

C_b2t = Rot_Mat_Fnc(z_h(7:9))';

F23 = -Skew_symmetric_matrix(C_b2t*u(1:3));   % See eq. (12.49) in Groves (2008) (however, here we use the tangent-frame, and not the earth-frame).

O=zeros(3);
I=eye(3);

C_e2t = Ce2t('N',init_lla(1),'E',init_lla(2));

init_t_frame = Ce2t('N',init_lla(1),'E',init_lla(2))*g2r('N',init_lla(1),'E',init_lla(2),init_lla(3));      % Calculate the "original" values of the origin in the t-frame.
[~,current_lat,~,~,~] = r2g(Ce2t('N',init_lla(1),'E',init_lla(2))'*(init_t_frame+z_h(1:3)));                % Calculate the current position in lla. 

g0 = norm(gravity(current_lat));                                                        % Calculate position dependent gravity. 
R0 = 6378137;                                                                           % See page 38 in Groves (2008).
e = 0.0818191908425;                                                                    % See page 38 in Groves (2008).
RE = R0/sqrt(1-e^2*sin(current_lat*180/pi)^2);                                          % See eq. (2.66) in Groves (2008).
r_eS2e = RE*sqrt(cos(current_lat*180/pi)^2+(1-e^2)*sin(current_lat*180/pi)^2);          % See eq. (2.89) in Groves (2008). 
F21 = 2*g0/(r_eS2e*norm(z_h(1:3))^2)*z_h(1:3)*z_h(1:3)';                                % See eq. (12.49) in Groves (2008) which shows the e-frame equivalent. 

% This is the matrix corresponding to eq. (12.48) in Groves (2008) (this is for the tangent-frame, and not the earth-frame).
Fc=[O I O O O O;
    F21 -2*C_e2t*simdata.Omega_ie2e*C_e2t' F23 C_b2t O O;
    O O -C_e2t*simdata.Omega_ie2e*C_e2t' O C_b2t O;
    O O O O O O;
    O O O O O O;
    O O O O O O];

% Calculate the process noise covariance matrix.
Q=zeros(12);
Q(1:3,1:3)=diag(simdata.sigma_acc).^2;
Q(4:6,4:6)=diag(simdata.sigma_gyro).^2;
Q(7:9,7:9)=simdata.sigma_acc_bias^2*eye(3);
Q(10:12,10:12)=simdata.sigma_gyro_bias^2*eye(3);
Q(13:15,13:15)=simdata.sigma_smph_process_noise^2*eye(3);

% Approximation of the discret time transition matrix
F=eye(18)+dT*Fc;

% Noise gain matrix
G=[O O O O O; C_b2t O O O O; O C_b2t O O O; O O I O O; O O O I O; O O O O I];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Correction of navigation state    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [z_out,C_b2s]=correct_navigation_estimates(z_in,dz,dPsi,C_b2s)

z_out = zeros(9,1);

z_out(1:6)=z_in(1:6)-dz(1:6);                                 % See eq. (5.93) in Groves (2008).

C_s2t = Rot_Mat_Fnc(z_in(7:9))';

dPsits2t = Skew_symmetric_matrix(dz(7:9));              
d_Cs2t = (eye(3)+dPsits2t);                              % See eq. (5.97) in Groves (2008).
C_s2t=d_Cs2t'*C_s2t;                                     % See eq. (5.95) in Groves (2008). Note that the lower row of eq. (5.96) is erroneous. 

dPsisb2s = Skew_symmetric_matrix(dPsi);              
d_Cb2s = (eye(3)+dPsisb2s);                              % See eq. (5.97) in Groves (2008).
C_b2s=d_Cb2s'*C_b2s;  

% Calculate yaw, pitch, roll. See eq. (2.17) in Groves (2008).
z_out(7)=atan2(C_s2t(3,2),C_s2t(3,3));                       % roll.
z_out(8)=-asin(C_s2t(3,1));                                  % pitch. 
z_out(9)=atan2(C_s2t(2,1),C_s2t(1,1));                       % yaw.

C_b2s = C_b2s*sqrtm(eye(3)/(C_b2s'*C_b2s)); % Find nearest orthogonal matrix.

z_out = real(z_out); 

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Gain calculation      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [K,H,R]=Gain(P,z_h,C_b2s)

global simdata;

% Measurement noise covariance matrix.
R=blkdiag(simdata.sigma_gps_pos_horiz^2*eye(2),simdata.sigma_gps_pos_vert^2,simdata.sigma_gps_speed^2,simdata.sigma_gps_speed^2/min(40,max((norm(z_h(4:6))^2),0.1)),simdata.sigma_vel_pseudo_obs_side^2,simdata.sigma_vel_pseudo_obs_down^2,simdata.pseudo_roll^2);

N1=length(P);
N2=length(R);

% Calculate observation matrix.
H=zeros(N2,N1);
H(1:3,1:3)=eye(3);

if(norm(z_h(4:5))>0)
    H(4,4) = z_h(4)/norm(z_h(4:5));
    H(4,5) = z_h(5)/norm(z_h(4:5));
else
    H(4,4) = 1/sqrt(2); 
    H(4,5) = 1/sqrt(2);
end
H(5,4) = -z_h(5)/norm(z_h(4:5))^2;
H(5,5) =  z_h(4)/norm(z_h(4:5))^2;

A = [zeros(2,1) eye(2)];
C_t2s = Rot_Mat_Fnc(z_h(7:9));

C_s2b = C_b2s';
% C_s2b = Rot_Mat_Fnc(rpy');

H(6:7,4:6) = A*C_s2b*C_t2s;
H(6:7,7:9) = A*C_s2b*C_t2s*Skew_symmetric_matrix(z_h(4:6));
H(6:7,16:18) = A*C_s2b*Skew_symmetric_matrix(C_t2s*z_h(4:6)); 
H(8,7) = 1;
H(8,16:18) = [1 0 0]*C_t2s';

% Calculate Kalman gain.
K=P*H'/(H*P*H'+R);

return

function  [K,H,R]=Gain2(P,z_h,C_b2s)

global simdata;

% Measurement noise covariance matrix.
R=blkdiag(simdata.sigma_vel_pseudo_obs_side^2,simdata.sigma_vel_pseudo_obs_down^2,simdata.pseudo_roll^2);

N1=length(P);
N2=length(R);

% Calculate observation matrix.
H=zeros(N2,N1);

A = [zeros(2,1) eye(2)];

C_t2s = Rot_Mat_Fnc(z_h(7:9));

% C_s2b = Rot_Mat_Fnc(rpy');
C_s2b = C_b2s';

H(1:2,4:6) = A*C_s2b*C_t2s;
H(1:2,7:9) = A*C_s2b*C_t2s*Skew_symmetric_matrix(z_h(4:6));
H(1:2,16:18) = A*C_s2b*Skew_symmetric_matrix(C_t2s*z_h(4:6)); 
H(3,7) = 1;
H(3,16:18) = [1 0 0]*C_t2s';

% Calculate Kalman gain.
K=P*H'/(H*P*H'+R);

atan2(C_b2s(3,2),C_b2s(3,3));

return
