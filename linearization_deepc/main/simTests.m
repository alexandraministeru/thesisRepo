% DeePC with wave forces preview, scaled signals and transfer functions.
%% Clean environment
clearvars;clc;close all;
rng('default')

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
addpath(genpath('functions'))

%% Get preview data
% Waves:
% % Hs=3, Tp=12
load('inputData\waveForces.mat');

% % % waves Hs=4.3, Tp=10
% outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs4p3_tp10_long\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
% [data, channels, units, headers] = ReadFASTtext(outFile);
% M_pitch = data(:,ismember(channels,'B1WvsMyi'));

M_pitch = M_pitch(1:20:end); % 1/0.05
clear data;

% Turbulent wind:n
load('inputData\turbWind_16mps_long.mat')
v = turbWind(1:20:end);

%% Load linearization
load('inputData\linDataWave.mat');
 
% %% Get equilibrium values of linearization point
% outFileOP = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.out';
% [dataOP, channelsOP, ~, ~] = ReadFASTtext(outFileOP);
% inputChannels = {'BldPitch1'}; % Collective blade pitch (deg)
% outputChannels = {'RotSpeed','PtfmPitch'};
% 
% % Input matrix
% nIn = length(inputChannels);
% U_OP = zeros(nIn,size(dataOP,1));
% 
% for idxCh = 1:length(inputChannels)
%     U_OP(idxCh,:) = dataOP(1:end,ismember(channelsOP,inputChannels{idxCh}))';
% end
% 
% timeSamples = size(dataOP,1);
% timeStep = dataOP(end,1)/(timeSamples - 1);
% timeWindow = 60; % seconds
% ssWindowIdx = timeWindow/timeStep;
% 
% u_opVal = mean(U_OP(:,end-ssWindowIdx:end),2);
% 
% % Find steady state operating value
% nChannels = length(outputChannels);
% y_opVal = zeros(nChannels,1);
% for idxCh=1:nChannels
%     y_opVal(idxCh) = getSSMean(dataOP, ssWindowIdx, channelsOP, outputChannels{idxCh});
% end

%% Find index of blade pitch, gen speed and wind speed
inputChannelsList = MBC.DescCntrlInpt;
outputChannelsList = MBC.DescOutput;

inputChannels = {'ED Extended input: collective blade-pitch command, rad', ...
    'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s',...
    'ED Platform Y moment, node 1, Nm'};

outputChannels = {'ED RotSpeed, (rpm)','ED PtfmPitch, (deg)','SrvD GenPwr, (kW)','ED TwrBsMyt, (kN-m)','ED RotTorq, (kN-m)','ED RootMyc1, (kN-m)'};

LTIsys_reduced = getReducedSS(MBC,LTIsys,inputChannels,outputChannels);

% Set input, output and state names
LTIsys_reduced.InputName = {'Collective blade pitch (\theta_c)', ...
    'Horizontal wind speed (v_0)', ...
    'Wave pitch moment (M_{pitch})'}; 
LTIsys_reduced.InputUnit = {'rad','m/s','Nm'};

LTIsys_reduced.OutputName = {'Rotor speed (\Omega)','Platform pitch angle (\Theta_P)','Generator power (P)','Tower base load (M_{y,T})', 'LSS load (M_{LSS})','Out-of-plane blade root load (M_{oop,B})'};
LTIsys_reduced.OutputUnit = {'rpm','deg','kW','kN-m','kN-m','kN-m'};

G_hat = LTIsys_reduced;

%% Scale transfer functions(G_hat, Gd_hat - unscaled; G, Gd scaled)
% Scaling factors
uhat_max = 5*(pi/180); % Maximum expected input (rad)
v0hat_max = 0.1*16; % Maximum expected wind disturbance (m/s)
MpitchHat_max = max(abs(M_pitch)); % Maximum expected wave pitch moment (Nm)

omegaHat_max = 0.15* 12.1; % Maximum expected generator speed error (10% around linearization OP) (rpm)
pltfmPitchHat_max = 1.5; % Maximum expected platform pitch angle (deg)

opVal_Pwr = 5000; %kW
opVal_Myt = 53318.25; %kN-m
opVal_LSS = 4179.43; %kN-m
opVal_Moop = 5743.38; %kN-m

pwrHat_max = 0.1*opVal_Pwr;
MytHat_max = 0.5*opVal_Myt;
LssHat_max = 0.5*opVal_LSS;
MoopHat_max = 0.5*opVal_Moop;

ehat_max = [omegaHat_max; pltfmPitchHat_max; pwrHat_max; MytHat_max; LssHat_max; MoopHat_max];

Du = diag([uhat_max v0hat_max MpitchHat_max]);
Dy = diag([omegaHat_max pltfmPitchHat_max pwrHat_max MytHat_max LssHat_max MoopHat_max]);

% Scaled transfer functions
G = Dy\G_hat*Du;

%% Discretize and get OL data
% Discretize system
Ts = 1;
G_d = c2d(G,Ts);

% Time vector
simTime = Ts*800;
t = 0:Ts:simTime-Ts;

% Generate PRBS input for persistency of excitation
u_hat_bladePitch = idinput(length(t),'PRBS',[0 1/10],[-2 2]); % in degrees
u_hat_bladePitch = u_hat_bladePitch*(pi/180); % in radians

% Wind input: setady wind
u_hat_windSpeed = zeros(size(u_hat_bladePitch));
% windData = 0.5.*randn(3000,1);
% u_hat_windSpeed = windData(1:size(u_hat_bladePitch,1));

% Wave disturbance input: pitch moment
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


controlParams.Std = Std;

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

dataOL.tsim = t;
dataOL.u = u_bladePitch;
dataOL.d = [u_windSpeed u_Mp];
dataOL.y = y;

%% DeePC parameters
N = 500; % lenght of data set
p = 20; % past data window
f = 30; % prediction window
Nbar = N-p-f+1;
i = 1;

controlParams.N = N;
controlParams.p = p;
controlParams.f = f;

%% Construct data matrices
% % Use preview information
previewAns = input('Use preview information? y/n[y]: ','s');

while not(isempty(previewAns)) && not(strcmp(previewAns,'y')) && ...
        not(strcmp(previewAns,'n'))
    disp('Invalid input.')
    previewAns = input('Use preview information? y/n[y]: ','s');
end

if isempty(previewAns) || strcmp(previewAns,'y')
    previewFlag = 1;
elseif strcmp(previewAns,'n')
    previewFlag = 0;
end

% previewFlag = 0;
y = y(:,1:2);

% Past data
data.Up = constructHankelMat(u_bladePitch,i,p,Nbar);
data.Yp = constructHankelMat(y,i,p,Nbar);

% Future data
data.Uf = constructHankelMat(u_bladePitch,i+p,f,Nbar);
data.Yf = constructHankelMat(y,i+p,f,Nbar);

disturbMat = u_Mp;
% disturbMat = [u_windSpeed u_Mp];

if previewFlag == 1
    data.Wp = constructHankelMat(disturbMat,i,p,Nbar); % past data
    data.Wf = constructHankelMat(disturbMat,i+p,f,Nbar); % future data
else
    data.Wp = [];
    data.Wf = [];
end

%% Set up control loop
kFinal = 800; % simulation steps
tsim = 0:Ts:Ts*(kFinal-1);
nInputs = size(data.Up,1)/p;
nOutputs = size(data.Yp,1)/p;
nDist = size(data.Wp,1)/p;

% % Generate reference trajectory
% ref = zeros(kFinal+f,1);
% % ref(200:end) = 100; % step in reference

% %MIMO reference:

ref_y1 = zeros(kFinal+f,1)';
% ref_y1(201:400) = 0.2/omegaHat_max*ones(200,1);
% ref_y1(601:800) = 0.2/omegaHat_max*ones(200,1);
% ref_y1 = zeros(kFinal+f,1)';
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
out = zeros(size(G_d.C,1),kFinal);

%% CL disturbances
% % % TurbWind
% load('inputData\turbWind_16mps_long.mat')
% v = turbWind(1:20:end);

% load('inputData\turbWind_16mps.mat') %turbulent wind obtained from a previous FAST simulation
% v = windData;

% % EOG
% load('inputData\eog_16mps.mat','Wind1VelX','Time')
% v = Wind1VelX(1:125:end);
% % v = interp1(Time,Wind1VelX,tsim)'; % resample with current sampling period
% v = v-16; % center around linearization point
% % v = [v;zeros(f,1)];
% % v = v-16;
% v = [zeros(100,1); v; zeros(kFinal,1)];

% Center and scale wind disturbance
v = v-16; % center around linearization point
v = v./v0hat_max;

v = zeros(kFinal+f,1)/v0hat_max; % Steady wind
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

weightOutputs = diag([100 0]); % no reference tracking for platform pitch angle, just constraints down below
% weightOutputs = 1e2*diag(1);n
controlParams.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
weightInputs= diag(1);
controlParams.R = kron(eye(f),weightInputs);

% Choose input bounds
controlParams.lbu = -5*(pi/180);
controlParams.lbu = controlParams.lbu/uhat_max;
controlParams.ubu = 5*(pi/180);
controlParams.ubu = controlParams.ubu/uhat_max;

% Input rate constraint
duDeg = 8; % deg/s
duRad = duDeg*(pi/180); % rad/sn
duRad = duRad/uhat_max;
controlParams.duf = duRad*Ts;

scalingOut = diag(Dy);
controlledOutputsScaling = scalingOut(1:2);

% Output bounds

controlParams.lby = [-1.8; -0.35]; % rpm, deg
controlParams.lby = controlParams.lby./controlledOutputsScaling;
controlParams.uby = [1.8; 0.35]; % rpm, deg
controlParams.uby = controlParams.uby./controlledOutputsScaling;

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

controlParams.previewNoiseStd = 0;

tic
for k=1:kFinal
    disp('Iteration: ')
    disp(k)

    % Reference trajectory
    rf = ref((k-1)*nOutputs+1:(k-1)*nOutputs+nOutputs*f);

    % Wind preview
    if previewFlag == 1
        % data.wf = v(k:k+f-1);
        data.wf = Mp(k:k+f-1) + controlParams.previewNoiseStd*randn(1);
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
    data.yini = [data.yini(nOutputs+1:end); out(1:2,k)];
    if previewFlag == 1
        data.wini = [data.wini(nDist+1:end); Mp(k)];
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
plot(tsim,ref_y1(1:kFinal)*omegaHat_max) % reference
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
% % RMSE and standard deviation
% startUpDelay = controlParams.f + 10;
% info.rmse = rmse((out(1,startUpDelay:end).*omegaHat_max)',ref_y1(startUpDelay:kFinal)');fprintf('\nRMSE: %d',info.rmse);
% info.std = std(out(1,startUpDelay:end).*omegaHat_max);fprintf('\nStandard deviation: %d',info.std);

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
info.normE = norm(out(1,50:end));


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

saveFlag = 1;

if saveFlag == 1
    simType = 'steadyWind_irrWaves_Hs3Tp12';
    % simType = 'steadyWind_irrWaves_Hs4p3Tp10';
    % simType = 'trubWind_irrWaves_Hs3Tp12';
    % simType = 'trubWind_irrWaves_Hs4p3Tp10';
    inName = '\theta_c (in deg)';
    % outName = 'Generator speed (in rpm)';
    outName = LTIsys_reduced.OutputName;
    scaledFlag = 1;
    tuning = 0;
    tunedVar = 'var';
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
    while isfile(['focusedDeePC\' fileName ext])
        if isstrprop(fileName(end),'digit')
            currIdx = str2double(fileName(end));
            fileName = [fileName(1:end-1) num2str(currIdx+1)];
        else
            fileName = [fileName num2str(1)]; %#ok<AGROW>
        end
    end

    % Save variables
    if scaledFlag == 1
        save(['focusedDeePC\' fileName ext],'inName','outName','tsim','kFinal','ref','Ts','controlParams', 'out', ...
            'uSeq','description','scaledFlag','Du','Dy','dataOL','info');
    else
        save(['focusedDeePC\' fileName ext],'inName','outName','tsim','kFinal','ref','Ts','controlParams', 'out', ...
            'uSeq','description','scaledFlag','info');
    end
end