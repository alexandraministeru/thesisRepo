%% Plotting
load('outputData\turb_np_noss.mat')
out1 = rotSpeed;
uSeq1 = uSeq;

load('outputData\turb_wpStatic_noss_1.mat')
out2 = rotSpeed;
uSeq2 = uSeq;

load('outputData\turb_wpDynamic_noss_1.mat')
out3 = rotSpeed;
uSeq3 = uSeq;
% %%%
% load('outputData\wind_static.mat')
% Vhub1 = V_hub;
% 
% load('outputData\wind_dynamic.mat')
% Vhub2 = V_hub;
% 
% figure
% plot(Vhub1(:,3))
% hold on
% plot(Vhub2(:,3))
% %%%

stepIdxs = find(ref);

figure
% subplot(2,1,1) 
plot(tsim,out1)
hold on
plot(tsim,out2)
hold on
plot(tsim,out3)
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
title('DeePC for reference tracking')
grid on
hold on
plot(tsim,ref(1:length(out1))) % reference
xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
legend('No preview','Pre-simulated preview','QBlade preview','Reference','Location','SouthEast')
set(gcf,'Color','White')
% subplot(2,1,2)
% plot(tsim,v(1:length(tsim)))


% Control input sequence
figure
plot(tsim,uSeq1)
hold on
plot(tsim,uSeq2)
hold on
plot(tsim,uSeq3)
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-5 25])
yline(0,'r--','LineWidth',1)
yline(22,'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input')
legend('No preview','Pre-simulated preview','QBlade preview','Location','SouthEast')
grid on
set(gcf,'Color','White')


