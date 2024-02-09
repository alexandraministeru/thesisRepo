% DeePC with wave forces preview, scaled signals and transfer functions.
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
% outFile = 'D:\Master\TUD\Y2\Thesis\matlab\fromAmr\FF_OpenFAST\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
% [data, channels, units, headers] = ReadFASTtext(outFile);
% plotChannels = {'RotSpeed','GenSpeed','BldPitch1','GenTq','GenPwr','Wave1Elev','B1WvsFxi','B1WvsMyi','PtfmPitch'};
% PlotFASToutput(outFile,[],[],plotChannels,1)
%
% F_surge = data(:,ismember(channels,'B1WvsFxi'));
% M_pitch = data(:,ismember(channels,'B1WvsMyi'));
% save('inputData\waveForces.mat','F_surge','M_pitch');

load('inputData\waveForces.mat');
% M_pitch = M_pitch(1:20:end); % 1/0.05

%% Load linearization
% FilePath = "..\5MW_OC3Spar_DLL_WTurb_WavesIrr";
% MtlbTlbxPath = "D:\Master\TUD\Y2\Thesis\matlab\repo\linearization_deepc\fromAmr\matlab-toolbox-main";
% [LTIsys, MBC, matData, FAST_linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath);
% cd ..\main
% save('inputData\linDataWave.mat','VTK','FAST_linData','LTIsys','matData','MBC');

% load('inputData\linDataWave.mat');
%% Find index of blade pitch, gen speed and wind speed
% inputChannelsList = MBC.DescCntrlInpt;
% outputChannelsList = MBC.DescOutput;
% 
% inputChannels = {'ED Extended input: collective blade-pitch command, rad', ...
%     'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s',...
%     'ED Platform Y moment, node 1, Nm'};
%     % 'ED Platform X force, node 1, N', ...
% 
% % outputChannels = {'ED GenSpeed, (rpm)'};
% outputChannels = {'ED RotSpeed, (rpm)',
%     'ED PtfmPitch, (deg)',
%     'SrvD GenPwr, (kW)'};
% 
% LTIsys_reduced = getReducedSS(MBC,LTIsys,inputChannels,outputChannels);
% 
% % Set input, output and state names
% LTIsys_reduced.InputName = {'Collective blade pitch (\theta_c)', ...
%     'Horizontal wind speed (v_0)', ...
%     'Wave pitch moment (M_{pitch})'}; 
% % 'Wave surge force (F_{surge})',
% 
% LTIsys_reduced.InputUnit = {'rad','m/s','Nm'};
% 
% % LTIsys_reduced.OutputName = {'Generator speed (\Omega_g)'};
% LTIsys_reduced.OutputName = {'Rotor speed (\Omega)','Platform pitch angle','Generator power'};
% 
% LTIsys_reduced.OutputUnit = {'rpm','deg','kW'};
% % 
% LTIsys_reduced.StateName = {'Platform horizontal surge translation', ...
% 'Platform pitch tilt rotation', ...
% '1st tower fore-aft bending', ...
% 'First time derivative of Platform horizontal surge translation', ...
% 'First time derivative of Platform pitch tilt rotation', ...
% 'First time derivative of 1st tower fore-aft bending', ...
% 'First time derivative of Variable speed generator', ...
% 'ExctnPtfmSg1' , ...
% 'ExctnPtfmSg2' , ...
% 'ExctnPtfmSg3' , ...
% 'ExctnPtfmSg4' , ...
% 'ExctnPtfmSg5' , ...
% 'ExctnPtfmSg6' , ...
% 'ExctnPtfmSg7' , ...
% 'ExctnPtfmSg8' , ...
% 'ExctnPtfmSg9' , ...
% 'ExctnPtfmSg10', ...
% 'ExctnPtfmSg11', ...
% 'ExctnPtfmSg12', ...
% 'ExctnPtfmSg13', ...
% 'ExctnPtfmSg14', ...
% 'ExctnPtfmP1'  , ...
% 'ExctnPtfmP2'  , ...
% 'ExctnPtfmP3'  , ...
% 'ExctnPtfmP4'  , ...
% 'ExctnPtfmP5'  , ...
% 'ExctnPtfmP6'  , ...
% 'ExctnPtfmP7'  , ...
% 'ExctnPtfmP8'  , ...
% 'RdtnPtfmSg1'  , ...
% 'RdtnPtfmSg2'  , ...
% 'RdtnPtfmSg3'  , ...
% 'RdtnPtfmSg4'  , ...
% 'RdtnPtfmP1'   , ...
% 'RdtnPtfmP2'   , ...
% 'RdtnPtfmP3'   , ...
% 'RdtnPtfmP4'};
% 
% LTIsys_reduced.StateUnit = {'m', 'rad', 'm', 'm/s', 'rad/s', 'm/s', 'rad/s',...
%     '','','','','','','','','','','','','','','','','','','','','','','', ...
%     '','','','','','',''};
% 
% LTIsys_reduced.Name = 'NREL 5MW linearization around 16mps';
% 
% % save('inputData\waveLTIsys_reduced_rotSpeed.mat','LTIsys_reduced');
% % save('inputData\waveLTIsys_reduced_rotSpeed_ptfmPitch.mat','LTIsys_reduced');
% save('inputData\waveLTIsys_reduced_rotSpeed_ptfmPitch_genPwr.mat','LTIsys_reduced');

% load('inputData\waveLTIsys_reduced_rotSpeed.mat')
load('inputData\waveLTIsys_reduced_rotSpeed_ptfmPitch.mat');

G_hat = LTIsys_reduced;

%% Scale transfer functions(G_hat, Gd_hat - unscaled; G, Gd scaled)
% Scaling factors
uhat_max = 5*(pi/180); % Maximum expected input (rad)
v0hat_max = 0.1*16; % Maximum expected wind disturbance (m/s)
MpitchHat_max = max(abs(M_pitch)); % Maximum expected wave pitch moment (Nm)

omegaHat_max = 0.15* 12.1; % Maximum expected generator speed error (10% around linearization OP) (rpm)
pltfmPitchHat_max = 1.5; % Maximum expected platform pitch angle (deg)

ehat_max = [omegaHat_max; pltfmPitchHat_max];

Du = diag([uhat_max v0hat_max MpitchHat_max]);
Dy = diag([omegaHat_max pltfmPitchHat_max]);

% Scaled transfer functions
G = Dy\G_hat*Du;

%% Discretize and get OL data
% Discretize system
Ts = 0.05;
G_d = c2d(G,Ts);

% Time vector
simTime = Ts*800;
t = 0:Ts:simTime;

% Generate PRBS input for persistency of excitation
u_hat_bladePitch = idinput(length(t),'PRBS',[0 1/10],[-2 2]); % in degrees
u_hat_bladePitch = u_hat_bladePitch*(pi/180); % in radians

% Wind input: setady wind
% u_hat_windSpeed = zeros(size(u_hat_bladePitch));
windData = 0.5.*randn(3000,1);
u_hat_windSpeed = windData(1:size(u_hat_bladePitch,1));

% Wave disturbance input: pitch moment
% u_hat_Mp = std(M_pitch).*randn(length(t),1);
u_hat_Mp = M_pitch(1:length(u_hat_windSpeed));

% Scale inputs
u_bladePitch = u_hat_bladePitch./uhat_max;
u_windSpeed = u_hat_windSpeed./v0hat_max;
u_Mp = u_hat_Mp./MpitchHat_max;

% Simulate system and get open-loop I/O data set
y = lsim(G_d,[u_bladePitch u_windSpeed u_Mp],t); % zero initial condition, around linearization point

% % Add measurement noise on output
% noiseAns = input('Add measurement noise? y/n[y]: ','s');
% 
% while not(isempty(noiseAns)) && not(strcmp(noiseAns,'y')) && ...
%         not(strcmp(noiseAns,'n'))
%     disp('Invalid input.')
%     noiseAns = input('Add measurement noise? y/n[y]: ','s');
% end

noiseAns = 'y';

if isempty(noiseAns) || strcmp(noiseAns,'y')
    noiseFlag = 1;
elseif strcmp(noiseAns,'n')
    noiseFlag = 0;
end

Std = 1e-3; % measurement noise standard deviation
Std = Std*noiseFlag;
y = y + Std.*randn(size(y));

figure
plot(t,u_bladePitch)
grid on
xlabel('Time (in s)')
ylabel('Blade pitch angle, scaled (-)')
title('PRBS input')
set(gcf,'Color','White')

% figure
% plot(t,u_Mp)
% grid on
% xlabel('Time (in s)')
% ylabel('Wave pitch moment on platform (in N-m)')
% title('WvsMy')
% set(gcf,'Color','White')

figure
plot(t,y(:,1))
grid on
xlabel('Time (in s)')
ylabel('Rotor speed, scaled (-)')
% ylabel('Generator speed (in rpm)')
title('OL response to PE input - \omega_g around linearization point')
set(gcf,'Color','White')

figure
plot(t,y(:,2))
grid on
xlabel('Time (in s)')
ylabel('Platform pitch angle, scaled (-)')
title('OL response to PE input - \theta_P around linearization point')
set(gcf,'Color','White')

%% DeePC parameters
% N = 500; % lenght of data set
% p = 20; % past data window
% f = 20; % prediction windowN = 600; % lenght of data set
N = 500;
p = 40; % past data window
f = 20; % prediction window
Nbar = N-p-f+1;
i = 1;

controlParams.N = N;
controlParams.p = p;
controlParams.f = f;

%% Construct data matrices
% % Use preview information
% previewAns = input('Use preview information? y/n[y]: ','s');
% 
% while not(isempty(previewAns)) && not(strcmp(previewAns,'y')) && ...
%         not(strcmp(previewAns,'n'))
%     disp('Invalid input.')
%     previewAns = input('Use preview information? y/n[y]: ','s');
% end
% 
% if isempty(previewAns) || strcmp(previewAns,'y')
%     previewFlag = 1;
% elseif strcmp(previewAns,'n')
%     previewFlag = 0;
% end

previewFlag = 1;

% Past data
data.Up = constructHankelMat(u_bladePitch,i,p,Nbar);
data.Yp = constructHankelMat(y,i,p,Nbar);

% Future data
data.Uf = constructHankelMat(u_bladePitch,i+p,f,Nbar);
data.Yf = constructHankelMat(y,i+p,f,Nbar);

% disturbMat = u_Mp;
% disturbMat = [u_windSpeed u_Mp];
disturbMat = u_windSpeed;

if previewFlag == 1
    data.Wp = constructHankelMat(disturbMat,i,p,Nbar); % past data
    data.Wf = constructHankelMat(disturbMat,i+p,f,Nbar); % future data
else
    data.Wp = [];
    data.Wf = [];
end

%% Set up control loop
kFinal = 5000; % simulation steps
tsim = 0:Ts:Ts*(kFinal-1);
nInputs = size(data.Up,1)/p;
nOutputs = size(data.Yp,1)/p;
nDist = size(data.Wp,1)/p;

% % Generate reference trajectory
% ref = zeros(kFinal+f,1);
% % ref(200:end) = 100; % step in reference

% %MIMO reference:
ref_y1 = zeros(kFinal+f,1)';
ref_y2 = zeros(kFinal+f,1)';
ref = [ref_y1; ref_y2];
ref = reshape(ref,[],1);

% Keep track of states
nStates_G = size(G_d.A,1);
x_G = zeros(nStates_G,kFinal+1);

% Set initial condition
x0_G = zeros(nStates_G,1);
x_G(:,1) = x0_G;

% Keep track of input and output sequences
uSeq = zeros(nInputs,kFinal);
out = zeros(nOutputs,kFinal);

%% CL disturbances
% % % TurbWind
load('inputData\turbWind_16mps_long.mat')
% v = turbWind(1:20:end);
v = turbWind;

% load('inputData\turbWind_16mps.mat') %turbulent wind obtained from a previous FAST simulation
% v = windData;
v = v-16; % center around linearization point


% % EOG
% load('inputData\eog_16mps.mat','Wind1VelX','Time')
% v = Wind1VelX(1:125:end);
% % v = interp1(Time,Wind1VelX,tsim)'; % resample with current sampling period
% v = v-16; % center around linearization point
% % v = [v;zeros(f,1)];
% % v = v-16;
% v = [zeros(100,1); v; zeros(kFinal,1)];

% Scale wind disturbance
v = v./v0hat_max;

% v = zeros(kFinal+f,1)/v0hat_max; % Steady wind
% Mp = M_pitch./MpitchHat_max ;
Mp = M_pitch(length(u_Mp)+1:end)./MpitchHat_max;
% Mp = zeros(size(Mp));

%% Solve the constrained optimization problem
% Use instrumental variables
% ivAns = input('Use instrumental variables? y/n[y]: ','s');
% 
% while not(isempty(ivAns)) && not(strcmp(ivAns,'y')) && ...
%         not(strcmp(ivAns,'n'))
%     disp('Invalid input.')
%     ivAns = input('Use instrumental variables? y/n[y]: ','s');
% end
% 
% if isempty(ivAns) || strcmp(ivAns,'y')
%     ivFlag = 1;
% elseif strcmp(ivAns,'n')
%     ivFlag = 0;
% end
ivFlag = 1;

% Define control weights:
%
% weightOutputs diagonal matrix of size l-by-l, where l is the number of
% output channels and the n-th element on the diagonal represents the
% weight for the corresponding n-th output

weightOutputs = diag([1e2 0]); % no reference tracking for platform pitch angle, just constraints down below
% weightOutputs = 1e2*diag(1);
controlParams.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
weightInputs= 1*diag(1);
controlParams.R = kron(eye(f),weightInputs);

% Choose input bounds
controlParams.lbu = -5*(pi/180);
controlParams.lbu = controlParams.lbu/uhat_max;
controlParams.ubu = 5*(pi/180);
controlParams.ubu = controlParams.ubu/uhat_max;

% Input rate constraint
duDeg = 8; % deg/s
duRad = duDeg*(pi/180); % rad/s
duRad = duRad/uhat_max;
controlParams.duf = duRad*Ts;

% Output bounds
controlParams.lby = [-1.8; -1.5]; % rpm, degrees (10 for rotSpeed bc I don't bound it now)
controlParams.lby = controlParams.lby./diag(Dy);
controlParams.uby = [1.8; 1.5];
controlParams.uby = controlParams.uby./diag(Dy);

% Choose optimization method
% method = input(['Optimization method: 1 - QP, ' '2 - SDP, ' '3 - NLP: ']);
method = 1;

%% Control loop
% Past data for prediction
data.uini = constructHankelMat(u_bladePitch,i+N-p,p,1);
data.yini = constructHankelMat(y,i+N-p,p,1);

if previewFlag == 1
    data.wini = constructHankelMat(disturbMat,i+N-p,p,1);
else
    data.wini = [];
end

tic
for k=1:kFinal
    disp('Iteration: ')
    disp(k)

    % Reference trajectory
    rf = ref((k-1)*nOutputs+1:(k-1)*nOutputs+nOutputs*f);

    % Wind preview
    if previewFlag == 1
        data.wf = v(k:k+f-1);
        % data.wf = Mp(k:k+f-1);
        % previewData = [v(k:k+f-1) Mp(k:k+f-1)];
        % data.wf = reshape(previewData',[],1);
    else
        data.wf = [];
    end

    % DeePC optimal control input
    uStar = deepc(data,rf,controlParams,method,ivFlag,previewFlag);
    uSeq(:,k) = min(max(controlParams.lbu, uStar), controlParams.ubu); % Saturation

    u = [uSeq(:,k);
        v(k);
        Mp(k)];

    % Apply optimal input, simulate output
    x_G(:,k+1) = G_d.A*x_G(:,k) + G_d.B*u;
    out(:,k) = G_d.C*x_G(:,k) + G_d.D*u + Std.*randn(size(out(:,k)));

    % Update past data with most recent I/O data
    data.uini = [data.uini(nInputs+1:end); uSeq(:,k)];
    data.yini = [data.yini(nOutputs+1:end); out(:,k)];
    if previewFlag == 1
        data.wini = [data.wini(nDist+1:end); v(k)];
        % data.wini = [data.wini(nDist+1:end); v(k); Mp(k)];
    end
end
toc

%% Plotting
stepIdxs = find(ref);

% Controlled output
figure
plot(tsim,out(1,:).*omegaHat_max)
xlabel('Time (in s)')
ylabel([LTIsys_reduced.OutputName{1} ' (in ' LTIsys_reduced.OutputUnit{1} ')'])
title('DeePC for reference tracking')
grid on
hold on
plot(tsim,ref(1:kFinal)) % reference
xline(Ts*f,'k--','Future window size')
yline(controlParams.lby(1)*omegaHat_max,'r--','LineWidth',1)
yline(controlParams.uby(1)*omegaHat_max,'r--','LineWidth',1)
% xline(Ts*stepIdxs(1),'k--','Reference step')
legend('Controlled output','Reference','Location','SouthEast')
set(gcf,'Color','White')

% Controlled output
figure
plot(tsim,out(2,:).*pltfmPitchHat_max)
xlabel('Time (in s)')
ylabel([LTIsys_reduced.OutputName{2} ' (in ' LTIsys_reduced.OutputUnit{2} ')'])
title('DeePC for reference tracking')
grid on
hold on
% plot(tsim,ref(1:kFinal)) % reference
xline(Ts*f,'k--','Future window size')
yline(controlParams.lby(2)*pltfmPitchHat_max,'r--','LineWidth',1)
yline(controlParams.uby(2)*pltfmPitchHat_max,'r--','LineWidth',1)
% xline(Ts*stepIdxs(1),'k--','Reference step')
legend('Platform pitch angle','Location','SouthEast')
ylim([-5 5])
set(gcf,'Color','White')

% Control input sequence
figure
plot(tsim,rad2deg(uSeq.*uhat_max))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-15 15])
yline(rad2deg(controlParams.lbu*uhat_max),'r--','LineWidth',1)
yline(rad2deg(controlParams.ubu*uhat_max),'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input')
grid on
set(gcf,'Color','White')

% Control input rate
figure
plot(tsim(2:end), rad2deg((diff(uSeq).*uhat_max).*(1/Ts)))
xlabel('Time (in s)')
ylabel('\theta_c rate (in deg/s)')
yline(duDeg,'r--','LineWidth',1)
yline(-duDeg,'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
ylim([-15 15])
title('Control input rate')
grid on
set(gcf,'Color','White')

%% Analyse output
% RMSE and standard deviation
startUpDelay = controlParams.f + 10;
info.rmse = rmse((out(1,startUpDelay:end).*omegaHat_max)',ref_y1(startUpDelay:kFinal)');fprintf('\nRMSE: %d',info.rmse);
info.std = std(out(1,startUpDelay:end).*omegaHat_max);fprintf('\nStandard deviation: %d',info.std);

% % FFT and PSD of signals
% fWaveMin = 0.05; % Hz
% fWaveMax = 0.3; % Hz
% 
% % Rotor speed
% [f_rotSpeed,fft_rotSpeed,psd_rotSpeed] = getFFT(Ts,out(1,:).*omegaHat_max);
% 
% figure
% sp1 = subplot(2,1,1);
% plot(f_rotSpeed(find(f_rotSpeed==0):end),fft_rotSpeed(find(f_rotSpeed==0):end))
% xlabel('Frequency (in Hz)')
% ylabel('Amplitude')
% title('FFT')
% xlim([0 0.4])
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
% grid on
% sp2 = subplot(2,1,2);
% plot(f_rotSpeed(find(f_rotSpeed==0):end),psd_rotSpeed(find(f_rotSpeed==0):end))
% xlabel('Frequency (in Hz)')
% ylabel('Power/Frequency (rpm^2/Hz)')
% title('PSD')
% xlim([0 0.4])
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
% grid on
% linkaxes([sp1,sp2],'x');
% sgtitle('Rotor speed')
% set(gcf,'Color','White')
% 
% % Blade pitch angle
% [f_pitchControlIn,fft_pitchControlIn,psd_pitchControlIn] = getFFT(Ts,rad2deg(uSeq.*uhat_max));
% 
% figure
% sp1 = subplot(2,1,1);
% plot(f_pitchControlIn(find(f_pitchControlIn==0):end),fft_pitchControlIn(find(f_pitchControlIn==0):end))
% xlabel('Frequency (in Hz)')
% ylabel('Amplitude')
% title('FFT')
% xlim([0 0.4])
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
% grid on
% sp2 = subplot(2,1,2);
% plot(f_pitchControlIn(find(f_pitchControlIn==0):end),psd_pitchControlIn(find(f_pitchControlIn==0):end))
% xlabel('Frequency (in Hz)')
% ylabel('Power/Frequency (deg^2/Hz)')
% title('PSD')
% xlim([0 0.4])
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
% grid on
% linkaxes([sp1,sp2],'x');
% sgtitle('Blade pitch angle command')
% set(gcf,'Color','White')

% Pitch actuator use
pitchingRate = rad2deg((diff(uSeq).*uhat_max).*(1/Ts));
timeWindow = tsim;
n = numel(timeWindow);
ADC = 0;

for idxADC=1:n-1
    ADC = ADC + (abs(pitchingRate(idxADC))/duDeg)*Ts;
end

info.ADC = ADC/(timeWindow(end) - timeWindow(1));fprintf('\nADC: %d',ADC);


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
    simType = 'bp_rotsp_ptfmp';
    inName = '\theta_c (in deg)';
    % outName = 'Generator speed (in rpm)';
    outName = LTIsys_reduced.OutputName;
    scaledFlag = 1;
    tuning = 1;
    tunedVar = 'N';
    ext = '.mat';

    % Specify type of simulation
    fileName = simType;

    % Set description
    if previewFlag == 1
        if exist('disturbMat','var') == 0
            description = 'preview';
        elseif size(disturbMat,2) == 2
            description = 'Fsg_Mp_preview';
        elseif isequal(disturbMat,u_Mp)
             description = 'Mp_preview';
        else
             description = 'Fsg_preview';
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
    while isfile(['tuningRep\' fileName ext])
        if isstrprop(fileName(end),'digit')
            currIdx = str2double(fileName(end));
            fileName = [fileName(1:end-1) num2str(currIdx+1)];
        else
            fileName = [fileName num2str(1)]; %#ok<AGROW>
        end
    end

    % Save variables
    if scaledFlag == 1
        save(['tuningRep\' fileName ext],'inName','outName','tsim','kFinal','ref','Ts','controlParams', 'out', ...
            'uSeq','description','scaledFlag','uhat_max','ehat_max','omegaHat_max','info');
    else
        save(['tuningRep\' fileName ext],'inName','outName','tsim','kFinal','ref','Ts','controlParams', 'out', ...
            'uSeq','description','scaledFlag','info');
    end
end