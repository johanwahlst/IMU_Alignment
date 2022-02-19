% This function plots the angle errors.

function plot_angle_errors_tb2t(Psi_tb2t_errors,Cov_tb2t_tot)

step = 1/20; 
step_fac = 400;
start_symb = 200;

%%%%%%%%%%
% plot 1 %
%%%%%%%%%%

close all
figure('units','normalized','position',[.1 .1 0.3 .15]);

y1 = sqrt(mean(Psi_tb2t_errors(1,:,:).^2,3))*180/pi; 
y2 = squeeze(mean(sqrt(Cov_tb2t_tot(1,1,:,:)),4))*180/pi;

clf
p1 = plot(1000,1000,'-rs');
hold on
p2 = plot(1000,1000,'--bo');
plot(0:step:step*(size(Psi_tb2t_errors,2)-1),y1,'--b')
plot(0+start_symb*step:step_fac*step:step*(size(Psi_tb2t_errors,2)-1),y1(1+start_symb:step_fac:end),'bo')
plot(0:step:step*(size(Psi_tb2t_errors,2)-1),y2,'-r')
plot(0+start_symb*step:step_fac*step:step*(size(Psi_tb2t_errors,2)-1),y2(1+start_symb:step_fac:end),'rs')

axis([0 120 0 10])
xlabel('xlab1');
ylabel('ylab1');
title('tit1');

leg_handle = legend([p1 p2],'leg1phonephonephon1','leg2phonephonephon1');
set(leg_handle,'Position',[0.6 0.62 0.25 0.2])

%%%%%%%%%%
% plot 2 %
%%%%%%%%%%

figure('units','normalized','position',[.1 .1 0.3 .15]);

y1 = sqrt(mean(Psi_tb2t_errors(2,:,:).^2,3))*180/pi; 
y2 = squeeze(mean(sqrt(Cov_tb2t_tot(2,2,:,:)),4))*180/pi;

clf
p1 = plot(1000,1000,'-rs');
hold on
p2 = plot(1000,1000,'--bo');
plot(0:step:step*(size(Psi_tb2t_errors,2)-1),y1,'--b')
plot(0+start_symb*step:step_fac*step:step*(size(Psi_tb2t_errors,2)-1),y1(1+start_symb:step_fac:end),'bo')
plot(0:step:step*(size(Psi_tb2t_errors,2)-1),y2,'-r')
plot(0+start_symb*step:step_fac*step:step*(size(Psi_tb2t_errors,2)-1),y2(1+start_symb:step_fac:end),'rs')

axis([0 120 0 10])
xlabel('xlab2');
ylabel('ylab2');
title('tit2');

leg_handle = legend([p1 p2],'leg1phonephonephon1','leg2phonephonephon1');
set(leg_handle,'Position',[0.6 0.62 0.25 0.2])

%%%%%%%%%%
% plot 3 %
%%%%%%%%%%

figure('units','normalized','position',[.1 .1 0.3 .15]);

y1 = sqrt(mean(Psi_tb2t_errors(3,:,:).^2,3))*180/pi; 
y2 = squeeze(mean(sqrt(Cov_tb2t_tot(3,3,:,:)),4))*180/pi;

clf
p1 = plot(1000,1000,'-rs');
hold on
p2 = plot(1000,1000,'--bo');
plot(0:step:step*(size(Psi_tb2t_errors,2)-1),y1,'--b')
plot(0+start_symb*step:step_fac*step:step*(size(Psi_tb2t_errors,2)-1),y1(1+start_symb:step_fac:end),'bo')
plot(0:step:step*(size(Psi_tb2t_errors,2)-1),y2,'-r')
plot(0+start_symb*step:step_fac*step:step*(size(Psi_tb2t_errors,2)-1),y2(1+start_symb:step_fac:end),'rs')
axis([0 120 0 10])

xlabel('xlab3');
ylabel('ylab3');
title('tit3');

leg_handle = legend([p1 p2],'leg1phonephonephon1','leg2phonephonephon1');
set(leg_handle,'Position',[0.6 0.62 0.25 0.2])


end

