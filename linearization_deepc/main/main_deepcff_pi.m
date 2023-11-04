% PI+DeePC with wave forces preview, scaled signals and transfer functions.
%% Clean environment
clearvars;clc;close all;
rng('default')

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
addpath(genpath('functions'))
% Insert full path to casadi lib
addpath(genpath('D:\Program Files\MATLAB\R2023a\casadi'));
import casadi.*
% YALMIP and solvers
addpath(genpath('D:\Program Files\mosek\10.1\toolbox\r2017aom')) %MOSEK
addpath(genpath('D:\Program Files\MATLAB\R2023a\sedumi')) % sedumi
addpath(genpath('D:\Program Files\MATLAB\R2023a\sdpt3')) % sdpt3
addpath(genpath('D:\Program Files\MATLAB\R2023a\yalmip')) % yalmip

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
MpitchHat_max = max(M_pitch); % Maximum expected wave pitch moment (Nm)
ehat_max = 0.1* 12.1; % Maximum expected generator speed error (10% around linearization OP) (rpm)

Du = diag([uhat_max v0hat_max MpitchHat_max]);
Dy = ehat_max;

% Scaled transfer functions
G = Dy\G_hat*Du;

%% Discretize and get OL data
Ts_wind = 0.05;
Ts_waves = 1;
[dataWind,controlParamsWind] = setupDeePCwind(G,Ts_wind);

Ts_ratio = Ts_waves/Ts_wind; % newTs/oldTs, only works when newTs is a multiple of oldTs
M_pitch = M_pitch(1:Ts_ratio:end);
[dataWaves,controlParamsWaves] = setupDeePCwave(G,Ts_waves,M_pitch);

% Discretize system
G_d = c2d(G,Ts_wind);

% Measurement noise standard deviation
Std = 5e-5;

%% Set up control loop
tFinal = 75; % simualtion length in seconds

kFinalWind = tFinal/Ts_wind; % simulation steps
kFinalWaves = tFinal/Ts_waves;
tsimWind = 0:Ts_wind:tFinal-Ts_wind;
tsimWaves= 0:Ts_waves:tFinal-Ts_waves;

nInputs = 1; % Control inputs
nOutputs = 1; % Controlled outputs
nDist = 1; % Disturbance channels, both controllers have 1 disturbance chanel

% Generate reference trajectory
refWind = zeros(kFinalWind + controlParamsWind.f,1);
refWaves = zeros(kFinalWaves + controlParamsWaves.f,1);

% Keep track of states
nStates_G = size(G_d.A,1);
x_G = zeros(nStates_G,kFinalWind+1);

% Set initial condition
x0_G = zeros(nStates_G,1);
x_G(:,1) = x0_G;

% Keep track of input and output sequences
uSeqFB = zeros(nInputs,kFinalWind);
uSeqFF = zeros(nInputs,kFinalWind);
uTotal = zeros(nInputs,kFinalWind);
out = zeros(nOutputs,kFinalWind);

%% CL disturbances
% Waves:
Mp = M_pitch./MpitchHat_max ;

% % % Wind:
% % % Turbulent wind
% load('inputData\turbWind_16mps_long.mat')
% v = [turbWind; turbWind; turbWind];
% v = v-16; % center around linearization point
% v = v./v0hat_max;

% % EOG
load('inputData\eog_16mps.mat','Wind1VelX','Time')
v = interp1(Time,Wind1VelX,0:Ts_wind:Ts_wind*(600-1))'; % resample with current sampling period
v = v-16;
v = [zeros(500,1); v ; zeros(2000,1)];
v = v./v0hat_max;

% v = zeros(20000,1);
% Mp = zeros(5000,1);


%% Solve the constrained optimization problem
ivFlag = 1;
previewFlagWind = 1;
previewFlagWaves = 1;

% Choose optimization method
% method = input(['Optimization method: 1 - QP, ' '2 - SDP, ' '3 - NLP: ']);
method = 1;

%% Control loop
tic
idxWaves = 1;

measOutput = 0;
filteredOutput = 0;
prevPitCom = 0;
prevIntErr = 0;
intErr = 0;

for idxWind=1:kFinalWind
    disp('Iteration: ')
    disp(idxWind)

    if mod(idxWind-1,Ts_ratio) == 0 

        % Reference trajectory: zero reference, disturbance rejection
        rfWaves = refWaves((idxWaves-1)*nOutputs+1:(idxWaves-1)*nOutputs+nOutputs*controlParamsWaves.f);

        % Wave preview
        dataWaves.wf = Mp(idxWaves:idxWaves+controlParamsWaves.f-1);
        uStarFF = deepc(dataWaves,rfWaves,controlParamsWaves,method,ivFlag,previewFlagWaves);
        uSeqFF(:,idxWind) = uStarFF(1);        

        idxWaves = idxWaves + 1;
    else
        uSeqFF(:,idxWind) = uSeqFF(:,idxWind-1);
    end

    if idxWaves> controlParamsWaves.f + 5
    % if idxWind >=2
        measOutput = out(:,idxWind-1);
        prevPitCom = uSeqFB(:,idxWind-1);
        prevIntErr = intErr;
    % end
        % %%%%%%%%%%%%%%%%%%%%%%%%% LPF %%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % LPF_cutoff = 0.157; % rad/s
        % Alpha = exp( -Ts_wind * LPF_cutoff );
        % 
        % % Apply the Filter
        % filteredOutput = (1.0 - Alpha)*measOutput + Alpha * filteredOutput;
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    % [uSeqFB(:,idxWind),intErr] = getPIcom(refWind(idxWind),filteredOutput,prevPitCom,prevIntErr,Ts_wind);
    [uSeqFB(:,idxWind),intErr] = getPIcom(refWind(idxWind),measOutput,prevPitCom,prevIntErr,Ts_wind);

    end

    uTotal(:,idxWind) = uSeqFB(:,idxWind) + uSeqFF(:,idxWind);
    uTotal(:,idxWind) = min( max(uTotal(:,idxWind), deg2rad(-10)/uhat_max ),  deg2rad(10)/uhat_max );

    u = [uTotal(:,idxWind);
        v(idxWind)
        Mp(idxWaves)];

    % Apply optimal input, simulate output
    x_G(:,idxWind+1) = G_d.A*x_G(:,idxWind) + G_d.B*u;
    out(:,idxWind) = G_d.C*x_G(:,idxWind) + G_d.D*u + Std.*randn(size(out(:,idxWind)));    

    % if mod(idxWind-1,Ts_ratio) == 0
    %     dataWaves.uini = [dataWaves.uini(nInputs+1:end); uSeqFF(:,idxWind)];        
    %     dataWaves.yini = [dataWaves.yini(nOutputs+1:end); out(:,idxWind)];       
    %     if previewFlagWaves == 1
    %         dataWaves.wini = [dataWaves.wini(nDist+1:end); Mp(idxWaves)];
    %     end
    % end
end
toc

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
plot(tsimWind,refWind(1:kFinalWind)) % reference
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
yline(rad2deg(controlParamsWind.lbu*uhat_max),'r--','LineWidth',1)
yline(rad2deg(controlParamsWind.ubu*uhat_max),'r--','LineWidth',1)
% xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input FB')
grid on
set(gcf,'Color','White')

% Control input sequence FF
figure
plot(tsimWind,rad2deg(uSeqFF.*uhat_max))
hold on
plot(tsimWind,rad2deg(uSeqFF.*uhat_max),'kx')
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-15 15])
yline(rad2deg(controlParamsWaves.lbu*uhat_max),'r--','LineWidth',1)
yline(rad2deg(controlParamsWaves.ubu*uhat_max),'r--','LineWidth',1)
% xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input FF')
grid on
set(gcf,'Color','White')

% Control input sequence
figure
plot(tsimWind,rad2deg((uTotal).*uhat_max))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-15 15])
yline(rad2deg(controlParamsWind.lbu*uhat_max),'r--','LineWidth',1)
yline(rad2deg(controlParamsWind.ubu*uhat_max),'r--','LineWidth',1)
% xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
title('Total control input')
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

% Control input rate
figure
plot(tsimWind(1:end-1), rad2deg((diff(uSeqFF).*uhat_max).*(1/Ts_waves)))
xlabel('Time (in s)')
ylabel('\theta_c rate (in deg/s)')
yline(8,'r--','LineWidth',1)
yline(-8,'r--','LineWidth',1)
% xline(Ts*f,'k--','Future window size')
ylim([-15 15])
title('Control input rate FF')
grid on
set(gcf,'Color','White')

% Control input rate
figure
plot(tsimWind(1:end-1), rad2deg((diff(uTotal).*uhat_max).*(1/Ts_wind)))
xlabel('Time (in s)')
ylabel('\theta_c rate (in deg/s)')
yline(8,'r--','LineWidth',1)
yline(-8,'r--','LineWidth',1)
% xline(Ts*f,'k--','Future window size')
ylim([-15 15])
title('Control input rate uTot')
grid on
set(gcf,'Color','White')


%% Save data
% |   simType   |  description  | optMethod | scaled |     tuning     |
% |-------------|---------------|-----------|--------|----------------|
% |    wind     | noPreview     |    QP     |   []   |       []       |
% |    wave     | (wind)preview |    SDP    | scaled |tuning*tunedVar*|
% |             | Fsg_preview   |    NLP    |        |                |
% |             | Mp_preview    |           |        |                |
% |             | Fsg_Mp_preview|           |        |                |
% ---------------------------------------------------------------------

% saveAns = input('Save data? y/n[y]: ','s');
% if isempty(saveAns) || strcmp(saveAns,'y')
%     saveFlag = 1;
% elseif strcmp(saveAns,'n')
%     saveFlag = 0;
% end

saveFlag = 0;

if saveFlag == 1
    simType = 'waveRotSpeed';
    inName = '\theta_c (in deg)';
    % outName = 'Generator speed (in rpm)';
    outName = 'Rotor speed (in rpm)';
    scaledFlag = 1;
    tuning = 0;
    tunedVar = 'N';
    ext = '.mat';

    % Specify type of simulation
    fileName = simType;

    % Set description
    if previewFlagWind == 1
        if exist('disturbMat','var') == 0
            description = 'preview';
        elseif size(disturbMat,2) == 2
            description = 'Fsg_Mp_preview';
        elseif isequal(disturbMat,u_Fsg)
            description = 'Fsg_preview';
        else
            description = 'Mp_preview';
        end
        fileName = [fileName '_' description];
    else
        description = 'noPreview';
        fileName = [fileName '_noPreview'];
    end

    % Set optimization method
    switch method
        case 1
            optMethod = 'QP';
        case 2
            optMethod = 'NLP';
        case 3
            optMethod = 'SDP';
    end
    fileName = [fileName '_' optMethod];

    % Set if scaled
    if scaledFlag == 1
        fileName = [fileName '_scaled'];
    end

    % Set if tuned
    if tuning == 1
        fileName = [fileName '_tuning' tunedVar];
    end

    % Check if file name already exists (can go up to filename19.mat)
    while isfile(['outputData2c\' fileName ext])
        if isstrprop(fileName(end),'digit')
            currIdx = str2double(fileName(end));
            fileName = [fileName(1:end-1) num2str(currIdx+1)];
        else
            fileName = [fileName num2str(1)]; %#ok<AGROW>
        end
    end

    % Save variables
    if scaledFlag == 1
        save(['outputData2c\' fileName ext],'inName','outName','tsimWind','kFinal','refWind','Ts','controlParams', 'out', ...
            'uSeq','description','scaledFlag','uhat_max','ehat_max');
    else
        save(['outputData2c\' fileName ext],'inName','outName','tsimWind','kFinal','refWind','Ts','controlParams', 'out', ...
            'uSeq','description','scaledFlag');
    end
end