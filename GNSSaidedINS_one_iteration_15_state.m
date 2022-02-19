% This function performs one iteration of a standard GNSS-aided INS.
% 
% Inputs:
%           - y: GNSS measurements (pos, vel, bear).
%           - u: IMU measurements (acc + gyro).
%           - index_vector: Specifying the mapping between IMU and GNSS data.
%           - dT: Time between samples. 
%           - init_lla: Lat-Long-Alt in origin of tangent frame.
%           - P_in: State covariance matrix.
%           - z_h: Navigation solution.
%           - du: Estimated IMU bias.
%
% Outputs: 
%           - z_h: Navigation solution.
%           - u_h: IMU measurement with estimated bias subtracted.
%           - P_out: State covariance matrix.
%           - pred_err: Prediction error.
%           - du_out: Estimated IMU bias.
%           - eps: Normalized innovations.

function [z_h,u_h,P_out,pred_err,du_out,eps]=GNSSaidedINS_one_iteration_15_state(y,u,index_vector,dT,init_lla,P_in,z_in,du)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Run one iteration of filter         %%
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
   [K,H,R]=Gain(P_in,z_in);                                     % Calculation of gain matrix. 

   bearing_h = atan2(z_in(5),z_in(4));                          % Calculating estimated bearing. 

   if(abs(bearing_h-(y(5)-2*pi))<abs(bearing_h-y(5)))           % Adjusting the measured and estimated bearing to lie in the same interval. We assume that 0<meas.bear<2pi and -pi<est.bear<2pi.
       y(5) = y(5)-2*pi;
   elseif(abs(bearing_h-(y(5)+2*pi))<abs(bearing_h-y(5)))       % Adjusting the measured and estimated bearing to lie in the same interval. We assume that 0<meas.bear<2pi and -pi<est.bear<2pi.
       y(5) = y(5)+2*pi;
   end
   pred_err=y(:)-[z_in(1:3); norm(z_in(4:5)); bearing_h];       % Calculation of prediction error.
   dx=[zeros(9,1); du]+K*pred_err;                              % Error estimation: delta_x(k) = [delta_z(k) delta_u(k)].

   dz=dx(1:9);                                                  % Navigation errors.
   du_out=dx(10:end);                                           % Sensor bias.
    
   eps = (sqrtm(H*P_in*H'+R))\pred_err;
   P_out=(eye(15)-K*H)*P_in;                                    % Update covariance matrix P(k|k).
   z_h=correct_navigation_estimates(z_in,dz);                   % Error compensation.
else
   du_out = du;
   z_h = z_in;
   P_out = P_in;
   pred_err = zeros(5,1);
   eps = zeros(5,1);
end

% Compensate for sensor bias.
u_h = u-du_out;                                                  % See page 152 in Groves (2008).

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
Fc=[O I O O O;
    F21 -2*C_e2t*simdata.Omega_ie2e*C_e2t' F23 C_b2t O;
    O O -C_e2t*simdata.Omega_ie2e*C_e2t' O C_b2t;
    O O O O O;
    O O O O O];

% Calculate the process noise covariance matrix.
Q=zeros(12);
Q(1:3,1:3)=diag(simdata.sigma_acc).^2;
Q(4:6,4:6)=diag(simdata.sigma_gyro).^2;
Q(7:9,7:9)=simdata.sigma_acc_bias^2*eye(3);
Q(10:12,10:12)=simdata.sigma_gyro_bias^2*eye(3);

% Approximation of the discret
% time transition matrics
F=eye(15)+dT*Fc;

% Noise gain matrix
G=[O O O O; C_b2t O O O; O C_b2t O O; O O I O; O O O I];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Correction of navigation state    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z_out=correct_navigation_estimates(z_in,dz)

z_out = zeros(9,1);

z_out(1:6)=z_in(1:6)-dz(1:6);                                 % See eq. (5.93) in Groves (2008).

C_b2t = Rot_Mat_Fnc(z_in(7:9))';

dPsitb2t = Skew_symmetric_matrix(dz(7:9));              
d_Cb2t = (eye(3)+dPsitb2t);                              % See eq. (5.97) in Groves (2008).
C_b2t=d_Cb2t'*C_b2t;                                     % See eq. (5.95) in Groves (2008). Note that the lower row of eq. (5.96) is erroneous. 
C_b2t = real(C_b2t);

% Calculate yaw, pitch, roll. See eq. (2.17) in Groves (2008).
z_out(7)=atan2(C_b2t(3,2),C_b2t(3,3));                       % roll.
z_out(8)=-asin(C_b2t(3,1));                                  % pitch. 
z_out(9)=atan2(C_b2t(2,1),C_b2t(1,1));                       % yaw.

z_out = real(z_out);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Gain calculation      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [K,H,R]=Gain(P,z_h)

global simdata;

% Measurement noise covariance matrix.
R=blkdiag(simdata.sigma_gps_pos_horiz^2*eye(2),simdata.sigma_gps_pos_vert^2,simdata.sigma_gps_speed^2,simdata.sigma_gps_speed^2/min(40,max((norm(z_h(4:6))^2),0.1)));

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

H = -H;

% Calculate Kalman gain.
K=P*H'/(H*P*H'+R);

return

