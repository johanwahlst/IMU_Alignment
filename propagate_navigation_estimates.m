% This function calculates the current position, velocity and attitude 
% based upon the previous navigation state and the input from the IMU. 
% 
% Inputs:
%               - z_in: Navigation state.
%               - u: IMU measurements.
%               - init_lla: the initial position measurement in latitude, longitude, altitude.
%               - dT: the current time step.
%
% Outputs:
%               - z_out: Navigation state.
%               - C_b2t: Rotation matrix.

function [z_out]=propagate_navigation_estimates(z_in,u,init_lla,dT)

global simdata;

z_out=zeros(9,1);

init_t_frame = Ce2t('N',init_lla(1),'E',init_lla(2))*g2r('N',init_lla(1),'E',init_lla(2),init_lla(3));      % Calculate the "original" values of the origin in the t-frame.
[~,current_lat,~,~,~] = r2g(Ce2t('N',init_lla(1),'E',init_lla(2))'*(init_t_frame+z_in(1:3)));               % Calculate the current position in lla. 

g_t=[0 0 gravity(current_lat)]'; % Calculate position dependent gravity.

% Calculate yaw, pitch, roll. See eq. (2.17) in Groves (2008).
C_b2t = Rot_Mat_Fnc(z_in(7:9))';

C_e2t = Ce2t('N',init_lla(1),'E',init_lla(2));

Omega_ib2b = Skew_symmetric_matrix(u(4:6));
Omega_tb2b = real(-C_b2t'*C_e2t*simdata.Omega_ie2e*C_e2t'*C_b2t+Omega_ib2b);
C_b2t = C_b2t*(2*eye(3)+Omega_tb2b*dT)/(2*eye(3)-Omega_tb2b*dT); % See eq. (15) in Skog (2005); "A Low-cost GPS Aided Inertial Navigation System for Vehicle Applications".

C_b2t = real(C_b2t);

z_out(7)=atan2(C_b2t(3,2),C_b2t(3,3));
z_out(8)=-asin(C_b2t(3,1));
z_out(9)=atan2(C_b2t(2,1),C_b2t(1,1));
    
f_t=C_b2t*u(1:3);

z_out(1:3) = z_in(1:3)+z_in(4:6)*dT+(f_t+g_t-2*C_e2t*simdata.Omega_ie2e*z_in(4:6))*dT^2/2;                  % See eq. (5.31) in Groves(2008).
z_out(4:6) = z_in(4:6)+(f_t+g_t-2*C_e2t*simdata.Omega_ie2e*z_in(4:6))*dT;                                   % See eq. (5.29) in Groves(2008).

end



