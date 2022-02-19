% This function integrates GNSS and IMU data in a GNSS-aided INS.
% 
% Inputs:
%           - y: GNSS measurements (pos, vel, bear).
%           - u: IMU measurements (acc + gyro).
%           - index_vector: Specifying the mapping between IMU and GNSS data.
%           - time: Timestamps.
%           - init: Lat-Long-Alt in origin of tangent frame..
%
% Outputs: 
%           - z_h: Navigation solution.
%           - u_h: IMU measurement with estimated bias subtracted.
%           - var: State variance. 
%           - pred_err: Prediction errors.

function [z_h, du, var, pred_err]=GNSSaidedINS_15_state(y,u,index_vector,time,init_lla)

global simdata;

N=length(u);
M=length(y);

simdata.N=N;
simdata.init_heading = y(5,1); % Initialize heading from measurement.

[P]=initialize_state_covariance;

% Allocate vecors
[z_h,u_h,var,du,Id,pred_err]=initialize_vectors(N,M,P,u,y(:,1),time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Run filter                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=2;
for k=1:N-1
    % A standard measurement update is performed if there are GPS measurements.
    if index_vector(k+1)
       [K,H]=Gain(P,z_h(:,k));                                          % Calculation of gain matrix. 

       bearing_h = atan2(z_h(5,k),z_h(4,k));                            % Calculating estimated bearing. 

       if(abs(bearing_h-(y(5,m)-2*pi))<abs(bearing_h-y(5,m)))           % Adjusting the measured and estimated bearing to lie in the same interval. We assume that 0<meas.bear<2pi and -pi<est.bear<2pi.
           y(5,m) = y(5,m)-2*pi;
       elseif(abs(bearing_h-(y(5,m)+2*pi))<abs(bearing_h-y(5,m))) 
           y(5,m) = y(5,m)+2*pi;
       end
       pred_err_pseudo_obs = [zeros(2,1) eye(2)]*Rot_Mat_Fnc(z_h(7:9,k))*z_h(4:6,k);
       pred_err(:,m)=[z_h(1:3,k); norm(z_h(4:5,k)); bearing_h; pred_err_pseudo_obs]-[y(:,m); zeros(2,1)];   % Calculation of prediction error.       dx=[zeros(9,1); du]+K*pred_err(:,m);                             % Error estimation: delta_x(k) = [delta_z(k) delta_u(k)].
       dx=[zeros(9,1); du]+K*pred_err(:,m);  
       m=m+1;                                                           % Increase counter of used GPS measurements.

       dz=dx(1:9);                                                      % Navigation errors.
       du=dx(10:end);                                                   % Sensor bias.

       P=(Id-K*H)*P;                                                    % Update covariance matrix P(k|k).
       z_h(:,k)=correct_navigation_estimates(z_h(:,k),dz);            % Error compensation.
    % Else, a measurement update is performed only using pseudo observations.
    else
        [K,H]=Gain2(P,z_h(:,k));

        % Calculation of prediction error.
        pred_err_pseudo_obs=[zeros(2,1) eye(2)]*Rot_Mat_Fnc(z_h(7:9,k))*z_h(4:6,k);

        dx=[zeros(9,1); du]+K*pred_err_pseudo_obs;                       % Error estimation: delta_x(k) = [delta_z(k) delta_u(k)].

        dz=dx(1:9);                                                      % Navigation errors.
        du=dx(10:end);                                                   % Sensor bias.

        P=(Id-K*H)*P;                                                    % Update covariance matrix P(k|k).
        z_h(:,k)=correct_navigation_estimates(z_h(:,k),dz);            % Error compensation.
    end
    
    % Save the estimated variances. 
    var(:,k+1)=diag(P);

    % Compensate for sensor bias.
    u_h(:,k)=u(:,k)-du(1:6);                                            % See page 152 in Groves (2008).

    % Calculate time step.
    dT = time(k+1)-time(k);

    % Navigation equation update x(k+1)
    z_h(:,k+1)=propagate_navigation_estimates(z_h(:,k),u_h(:,k),init_lla,dT);

    % Calculate filter matrices.
    [F,G,Q]=calculate_filter_matrices(z_h(:,k+1),u_h(:,k),init_lla,dT);

    % Time update of estimated variance P(k+1|k), P=F*P*F'+G*Q*G';
    P=F*P*F'+G*Q*G'*dT; % See eq. (4.111) in Farrell (2008) for the approximation giving the second term.
end

u_h(:,end)=u(:,end)-du(1:6);    % Compensate for sensor bias at the last time step.

du = u-u_h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          SUB-FUNCTIONS                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Vector allocation       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [z_h,u_h,var,du,Id,pred_err]=initialize_vectors(N,M,P,u,y,time)

% Estimated state vectors.
z_h=zeros(9,N);
u_h=zeros(6,N);

z_h(1:3,1)=y(1:3,1); % Initialize position.
z_h(4:5,1) = y(4,1)*[cos(y(5,1)); sin(y(5,1))]; % Initialize velocity. 

% Vector of estimated variances.
var=zeros(15,N);

% Estimated accelerometer and gyroscope bias.
du=zeros(6,1);

% Prediction errors.
pred_err = zeros(7,M);

Id=eye(15);

% Calculate the initial roll and pitch. See eq. (5.89) in Groves (2008).
mean_force=mean(u(:,time<time(1)+1),2);
z_h(7,1)=atan2(-mean_force(2),-mean_force(3));
z_h(8,1)=atan2(mean_force(1),norm(mean_force(2:3)));
z_h(9,1)=y(5);                   

% Save the first variance.
var(:,1)=diag(P);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Initialize filter      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P]=initialize_state_covariance

global simdata;

% Covariance matrix.
P=zeros(15);

% General values for the initial covariance matrix P.
P(1:3,1:3)=simdata.init_sd(1)^2*eye(3);
P(4:6,4:6)=simdata.init_sd(2)^2*eye(3);
P(7:9,7:9)=diag(simdata.init_sd(3:5)).^2;
P(10:12,10:12)=simdata.init_sd(6)^2*eye(3);
P(13:15,13:15)=simdata.init_sd(7)^2*eye(3);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate_filter_matrices  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,G,Q]=calculate_filter_matrices(z_h,u,init_lla,dT)

global simdata

% Calculate rotation matrix.
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

% Calculate rotation matrix. 
C_b2t = Rot_Mat_Fnc(z_in(7:9))';

z_out(1:6)=z_in(1:6)-dz(1:6);                                % See eq. (5.93) in Groves (2008).

dPsitb2t = Skew_symmetric_matrix(dz(7:9));              
d_Cb2t = (eye(3)+dPsitb2t);                                  % See eq. (5.97) in Groves (2008).
C_b2t=d_Cb2t'*C_b2t;                                         % See eq. (5.95) in Groves (2008). Note that the lower row of eq. (5.96) is erroneous. 

% Calculate yaw, pitch, roll. See eq. (2.17) in Groves (2008).
z_out(7)=atan2(C_b2t(3,2),C_b2t(3,3));                       % roll.
z_out(8)=-asin(C_b2t(3,1));                                  % pitch. 
z_out(9)=atan2(C_b2t(2,1),C_b2t(1,1));                       % yaw.

z_out = real(z_out);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Gain calculation      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [K,H]=Gain(P,z_h)

global simdata;

% Measurement noise covariance matrix.
R=blkdiag(simdata.sigma_gps_pos_horiz^2*eye(2),simdata.sigma_gps_pos_vert^2,simdata.sigma_gps_speed^2,simdata.sigma_gps_speed^2/max((norm(z_h(4:6))^2),0.1),simdata.sigma_vel_pseudo_obs_side^2,simdata.sigma_vel_pseudo_obs_down^2);
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

if(norm(z_h(4:5)))
    H(5,4) = -z_h(5)/norm(z_h(4:5))^2;
    H(5,5) =  z_h(4)/norm(z_h(4:5))^2;
else
    H(5,4) = -1/2;
    H(5,5) = 1/2;
end

A = [zeros(2,1) eye(2)];
C_t2b = Rot_Mat_Fnc(z_h(7:9));

H(6:7,4:6) = A*C_t2b;
H(6:7,7:9) = A*C_t2b*Skew_symmetric_matrix(z_h(4:6));

% Calculate Kalman gain.
K=P*H'/(H*P*H'+R);

return

function  [K,H]=Gain2(P,z_h)

global simdata;

% Measurement noise covariance matrix.
R=blkdiag(simdata.sigma_vel_pseudo_obs_side^2,simdata.sigma_vel_pseudo_obs_down^2);

N1=length(P);
N2=length(R);

% Calculate observation matrix.
H=zeros(N2,N1);

A = [zeros(2,1) eye(2)];

C_t2b = Rot_Mat_Fnc(z_h(7:9));

H(1:2,4:6) = A*C_t2b;
H(1:2,7:9) = A*C_t2b*Skew_symmetric_matrix(z_h(4:6));

% Calculate Kalman gain.
K=P*H'/(H*P*H'+R);

return

