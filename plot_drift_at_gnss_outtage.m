% This function produces plots of the position drifts.

function plot_drift_at_gnss_outtage(pos_error1,pos_errors1,error_time1,pos_error2,pos_errors2,error_time2,pos_error3,pos_errors3,error_time3)

%%%%%%%%%%
% plot 1 %
%%%%%%%%%%

err11 = zeros(1,ceil(max(max(2*error_time1))));
err12 = zeros(1,ceil(max(max(2*error_time1))));
numb_of_meas = zeros(1,ceil(max(max(2*error_time1))));
for n = 1:size(pos_error1,1)
    for m = 1:size(pos_error1,2)
        ind = floor(error_time1(n,m)*2)+1;
        err11(ind) = err11(ind) + pos_error1(n,m);
        err12(ind) = err12(ind) + pos_errors1(n,m);
        numb_of_meas(ind) = numb_of_meas(ind) + 1; 
    end
end

err11 = err11./numb_of_meas; 
err12 = err12./numb_of_meas;

close all
figure('units','normalized','position',[.1 .1 0.3 .15]);

clf
p1 = plot(1000,1000,'--bo');
hold on
p2 = plot(1000,1000,'-rs');
plot(0:0.5:60,err11(1:121),'--b')
plot(0+2:5:60-2,err11(1+4:10:121-4),'ob')
plot(0:0.5:60,err12(1:121),'r')
plot(0+2:5:60-2,err12(1+4:10:121-4),'sr')
axis([0 60 0 800])
set(gca,'XTick',0:10:60);
set(gca,'YTick',0:200:800);
xlabel('xlabel1');
ylabel('ylabel1');
title('title1');
leg_handle = legend([p2 p1],'leg1phonephonephonephoneph1','leg2phonephonephonephoneph1');
set(leg_handle,'Position',[0.238 0.62 0.25 0.2])

%%%%%%%%%%
% plot 2 %
%%%%%%%%%%

err21 = zeros(1,ceil(max(max(2*error_time2))));
err22 = zeros(1,ceil(max(max(2*error_time2))));
numb_of_meas = zeros(1,ceil(max(max(2*error_time2))));
for n = 1:size(pos_error2,1)
    for m = 1:size(pos_error2,2)
        ind = floor(error_time2(n,m)*2)+1;
        err21(ind) = err21(ind) + pos_error2(n,m);
        err22(ind) = err22(ind) + pos_errors2(n,m);
        numb_of_meas(ind) = numb_of_meas(ind) + 1; 
    end
end

err21 = err21./numb_of_meas; 
err22 = err22./numb_of_meas;

figure('units','normalized','position',[.1 .1 0.3 .15]);

clf
p1 = plot(1000,1000,'--bo');
hold on
p2 = plot(1000,1000,'-rs');
plot(0:0.5:60,err21(1:121),'--b')
plot(0+2:5:60-2,err21(1+4:10:121-4),'ob')
plot(0:0.5:60,err22(1:121),'r')
plot(0+2:5:60-2,err22(1+4:10:121-4),'sr')
axis([0 60 0 800])
set(gca,'XTick',0:10:60);
set(gca,'YTick',0:200:800);
xlabel('xlabel2');
ylabel('ylabel2');
title('title2');
leg_handle = legend([p2 p1],'leg1phonephonephonephoneph2','leg2phonephonephonephoneph2');
set(leg_handle,'Position',[0.238 0.62 0.25 0.2])

%%%%%%%%%%
% plot 3 %
%%%%%%%%%%

err31 = zeros(1,ceil(max(max(2*error_time3))));
err32 = zeros(1,ceil(max(max(2*error_time3))));
numb_of_meas = zeros(1,ceil(max(max(2*error_time3))));
for n = 1:size(pos_error3,1)
    for m = 1:size(pos_error3,2)
        ind = floor(error_time3(n,m)*2)+1;
        err31(ind) = err31(ind) + pos_error3(n,m);
        err32(ind) = err32(ind) + pos_errors3(n,m);
        numb_of_meas(ind) = numb_of_meas(ind) + 1; 
    end
end

err31 = err31./numb_of_meas; 
err32 = err32./numb_of_meas;

figure('units','normalized','position',[.1 .1 0.3 .15]);

clf
p1 = plot(1000,1000,'--bo');
hold on
p2 = plot(1000,1000,'-rs');
plot(0:0.5:60,err31(1:121),'--b')
plot(0+2:5:60-2,err31(1+4:10:121-4),'ob')
plot(0:0.5:60,err32(1:121),'r')
plot(0+2:5:60-2,err32(1+4:10:121-4),'sr')
axis([0 60 0 800])
set(gca,'XTick',0:10:60);
set(gca,'YTick',0:200:800);
xlabel('xlabel3');
ylabel('ylabel3');
title('title3');
leg_handle = legend([p2 p1],'leg1phonephonephonephoneph3','leg2phonephonephonephoneph3');
set(leg_handle,'Position',[0.238 0.62 0.25 0.2])

end

