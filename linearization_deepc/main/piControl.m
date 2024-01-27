% PI+DeePC with wave forces preview, scaled signals and transfer functions.
%% Clean environment
clearvars;clc;close all;
rng('default')

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
addpath(genpath('functions'))

%% Get preview data
load('inputData\waveForces.mat','M_pitch');

%% Load linearization
load('inputData\linDataWave.mat');

%% Find index of blade pitch, gen speed and wind speed
load('inputData\waveLTIsys_reduced_rotSpeed.mat')
G_hat = LTIsys_reduced;

%% Scale transfer functions(G_hat, Gd_hat - unscaled; G, Gd scaled)
% Scaling factors
uhat_max = 10*(pi/180); % Maximum expected input (rad)
v0hat_max = 0.1*16; % Maximum expected wind disturbance (m/s)
MpitchHat_max = max(abs(M_pitch)); % Maximum expected wave pitch moment (Nm)
ehat_max = 0.1* 12.1; % Maximum expected generator speed error (10% around linearization OP) (rpm)

Du = diag([uhat_max v0hat_max MpitchHat_max]);
Dy = ehat_max;

% Scaled transfer functions
G = Dy\G_hat*Du;

%% Discretize and get OL data
Ts_wind = 0.05;

% Discretize system
G_d = c2d(G,Ts_wind);

% Measurement noise standard deviation
Std = 1e-5;

%% CL data collection
% Time vector
simTime = 120;
t = 0:Ts_wind:simTime-Ts_wind;
kFinal = simTime/0.05;

ref = zeros(length(t),1);

% % Turbulent wind
load('inputData\turbWind_16mps.mat')
v = [windData; windData; windData];
v = v-16; % center around linearization point
v = v./v0hat_max;

% Mp = zeros(kFinal,1);
Mp = M_pitch./MpitchHat_max;

uSeqFB = zeros(1,kFinal);
uTotal = zeros(1,kFinal);
out = zeros(1,kFinal);

% Keep track of states
nStates_G = size(G_d.A,1);
x_G = zeros(nStates_G,kFinal+1);

% Set initial condition
x0_G = zeros(nStates_G,1);
x_G(:,1) = x0_G;

measOutput = 0;
filteredOutput = 0;
prevPitCom = 0;
prevIntErr = 0;
intErr = 0;

% Closed loop
for idxCL = 1:kFinal
    disp('Iteration: ')
    disp(idxCL)
    if idxCL>1
        measOutput = out(:,idxCL-1);
        prevPitCom = uSeqFB(:,idxCL-1);
        prevIntErr = intErr;
    end

    if idxCL > 25
        filteredOutput = lowpass(out(:,idxCL-25:idxCL),3e-3,Ts_wind);
    else
        filteredOutput = lowpass(out(:,1:idxCL),3e-3,Ts_wind);
    end
    [uSeqFB(:,idxCL),intErr] = getPIcom(ref(idxCL),filteredOutput(end),prevPitCom,prevIntErr,Ts_wind);

    % Alpha= exp(-Ts_wind *0.25);
    % filteredOutput = ( 1.0 - Alpha )*measOutput + Alpha * filteredOutput;
    % [uSeqFB(:,idxCL),intErr] = getPIcom(ref(idxCL),filteredOutput,prevPitCom,prevIntErr,Ts_wind);    

    % [uSeqFB(:,idxCL),intErr] = getPIcom(ref(idxCL),measOutput,prevPitCom,prevIntErr,Ts_wind);

    uTotal(:,idxCL) = uSeqFB(:,idxCL);
    uTotal(:,idxCL) = min( max(uTotal(:,idxCL), deg2rad(-10)/uhat_max ),  deg2rad(10)/uhat_max );

    u = [uTotal(:,idxCL);
        v(idxCL);
        Mp(idxCL)];

    % Apply optimal input, simulate output
    x_G(:,idxCL+1) = G_d.A*x_G(:,idxCL) + G_d.B*u;
    out(:,idxCL) = G_d.C*x_G(:,idxCL) + G_d.D*u + Std.*randn(size(out(:,idxCL)));
end

tsimWind = t;



%% Plotting
% stepIdxs = find(refWind);

% Controlled output
figure
plot(tsimWind,out.*ehat_max)
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
% ylabel('Generator speed (in rpm)')
title('DeePC for reference tracking')
grid on
hold on
plot(tsimWind,ref(1:kFinal)) % reference
% xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
legend('Controlled output','Reference','Location','SouthEast')
set(gcf,'Color','White')

% Control input sequence FB
figure
plot(tsimWind,rad2deg(uSeqFB.*uhat_max))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-15 15])
yline(-10,'r--','LineWidth',1)
yline(10,'r--','LineWidth',1)
% xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input FB')
grid on
set(gcf,'Color','White')


% Control input rate
figure
plot(tsimWind(1:end-1), rad2deg((diff(uSeqFB).*uhat_max).*(1/Ts_wind)))
xlabel('Time (in s)')
ylabel('\theta_c rate (in deg/s)')
yline(8,'r--','LineWidth',1)
yline(-8,'r--','LineWidth',1)
% xline(Ts*f,'k--','Future window size')
ylim([-15 15])
title('Control input rate FB')
grid on
set(gcf,'Color','White')



