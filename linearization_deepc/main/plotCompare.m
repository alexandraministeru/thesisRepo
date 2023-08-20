%% Plotting
load('outputData\uRateConstr_quadprog_np.mat')
out1 = out;
uSeq1 = uSeq;

load('outputData\uRateConstr_quadprogIP.mat')
out2 = out;
uSeq2 = uSeq;

stepIdxs = find(ref);

figure
plot(tsim,out1)
hold on
plot(tsim,out2)
xlabel('Time (in s)')
ylabel('Generator speed (in rpm)')
title('DeePC for reference tracking')
grid on
hold on
plot(tsim,ref(1:kFinal)) % reference
xline(Ts*f,'k--','Future window size')
xline(Ts*stepIdxs(1),'k--','Reference step')
legend('No preview','With preview','Reference','Location','SouthEast')
set(gcf,'Color','White')

% Control input sequence
figure
plot(tsim,uSeq1*(180/pi))
hold on
plot(tsim,uSeq2*(180/pi))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-15 15])
yline(-10,'r--','LineWidth',1)
yline(10,'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input')
legend('No preview','With preview','Reference','Location','SouthEast')
grid on
set(gcf,'Color','White')


