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
Ts_waves = 1;
Ts_ratio = Ts_waves/Ts_wind; % newTs/oldTs, only works when newTs is a multiple of oldTs

% Discretize system
G_d = c2d(G,Ts_wind);

% Measurement noise standard deviation
Std = 1e-5;

%% CL data collection
% Time vector
simTime = 500;
tWind = 0:Ts_wind:simTime-Ts_wind;
tWaves = 0:Ts_waves:simTime-Ts_waves;

% Generate PRBS input for persistency of excitation
ref_hat = zeros(length(tWind),1); % rpm
windData = 0.5.*randn(length(tWind),1);
u_hat_windSpeed = windData(1:length(tWind));
u_hat_Mp = M_pitch(1:length(tWind));

u_hat_bladePitch = idinput(length(tWaves),'PRBS',[0 1/10],[-1 1]); % in degrees
u_hat_bladePitch = u_hat_bladePitch*(pi/180); % in radians

u_bladePitch = u_hat_bladePitch./uhat_max;
u_windSpeed = u_hat_windSpeed./v0hat_max;
u_Mp = u_hat_Mp./MpitchHat_max;
ref = ref_hat./ehat_max;

u_CL = zeros(length(tWaves),1);
u_CL_PI = zeros(length(tWaves),1);
d_CL = zeros(length(tWaves),1);
y_CL = zeros(length(tWaves),1);

uSeqFB = zeros(1,length(tWaves));
uTotal = zeros(1,length(tWind));
out = zeros(1,length(tWind));

% Keep track of states
nStates_G = size(G_d.A,1);
x_G = zeros(nStates_G,length(tWind)+1);

% Set initial condition
x0_G = zeros(nStates_G,1);
x_G(:,1) = x0_G;

measOutput = 0;
filteredOutput = 0;
prevPitCom = 0;
prevIntErr = 0;
intErr = 0;

% Closed loop
for idxCL = 1:length(tWind)
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

    idxWaves = mod(idxCL-1,20) + 1;
    uTotal(:,idxCL) = uSeqFB(:,idxCL) + u_bladePitch(idxWaves);
    uTotal(:,idxCL) = min( max(uTotal(:,idxCL), deg2rad(-10)/uhat_max ),  deg2rad(10)/uhat_max );

    u = [uTotal(:,idxCL);
        u_windSpeed(idxCL);
        u_Mp(idxCL)];

    % Apply optimal input, simulate output
    x_G(:,idxCL+1) = G_d.A*x_G(:,idxCL) + G_d.B*u;
    out(:,idxCL) = G_d.C*x_G(:,idxCL) + G_d.D*u + Std.*randn(size(out(:,idxCL)));

    
    if mod(idxCL-1,20) == 0
        idxWaves = floor((idxCL-1)/20)+1;
        u_CL(idxWaves) = u_bladePitch(idxWaves);
        u_CL_PI(idxWaves) = uSeqFB(:,idxCL);
        y_CL(idxWaves) = out(:,idxCL);
        d_CL(idxWaves) = u_Mp(idxCL);
    end

end

save('inputData\clPI.mat','u_CL',"d_CL","y_CL","u_CL_PI");

%% DeePC parameters
N = 500; % lenght of data set
p = 20; % past data window
f = 20; % prediction windowN = 600; % lenght of data set
Nbar = N-p-f+1;
i = 1;

controlParams.N = N;
controlParams.p = p;
controlParams.f = f;
Ts=1;

%% Construct data matrices
% Use preview information
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

uMat = [u_CL u_CL_PI];
% uMat = u_CL;

% Past data
data.Up = constructHankelMat(uMat,i,p,Nbar);
data.Yp = constructHankelMat(y_CL,i,p,Nbar);

% Future data
data.Uf = constructHankelMat(uMat,i+p,f,Nbar);
data.Yf = constructHankelMat(y_CL,i+p,f,Nbar);

% disturbMat = u_windSpeed;
disturbMat = d_CL;
% disturbMat = [u_windSpeed u_Mp];

if previewFlag == 1
    data.Wp = constructHankelMat(disturbMat,i,p,Nbar); % past data
    data.Wf = constructHankelMat(disturbMat,i+p,f,Nbar); % future data
else
    data.Wp = [];
    data.Wf = [];
end

%% Solve the constrained optimization problem
ivFlag = 1;

% Define control weights:
%
% weightOutputs diagonal matrix of size l-by-l, where l is the number of
% output channels and the n-th element on the diagonal represents the
% weight for the corresponding n-th output

weightOutputs = 1e2; % no reference tracking for platform pitch angle, just constraints down below
% weightOutputs = 1e2*diag(1);
controlParams.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
weightInputs= 1*diag([1 0]);
% weightInputs= 1*diag(1);
controlParams.R = kron(eye(f),weightInputs);

% Choose input bounds
controlParams.lbu = -10*(pi/180);
controlParams.lbu = controlParams.lbu/uhat_max;
controlParams.lbu = [controlParams.lbu;controlParams.lbu];
controlParams.ubu = 10*(pi/180);
controlParams.ubu = controlParams.ubu/uhat_max;
controlParams.ubu = [controlParams.ubu;controlParams.ubu];

% Input rate constraint
duDeg = 8; % deg/s
duRad = duDeg*(pi/180); % rad/s
duRad = duRad/uhat_max;
controlParams.duf = duRad*Ts;

controlParams.duf = [controlParams.duf; controlParams.duf];

% Choose optimization method
% method = input(['Optimization method: 1 - QP, ' '2 - SDP, ' '3 - NLP: ']);
method = 1;

%% Control loop
% Past data for prediction
data.uini = constructHankelMat(uMat,i+N-p,p,1);
data.yini = constructHankelMat(y_CL,i+N-p,p,1);
data.wini = constructHankelMat(d_CL,i+N-p,p,1);

%% Set up control loop
tFinal = 100; % simualtion length in seconds

kFinalWind = tFinal/Ts_wind; % simulation steps
kFinalWaves = tFinal/Ts_waves;
tsimWind = 0:Ts_wind:tFinal-Ts_wind;
tsimWaves= 0:Ts_waves:tFinal-Ts_waves;

nInputs = size(data.Up,1)/controlParams.p; % Control inputs
nOutputs = 1; % Controlled outputs
nDist = 1; % Disturbance channels, both controllers have 1 disturbance chanel

% Generate reference trajectory
refWind = zeros(kFinalWind + f,1);
refWaves = zeros(kFinalWaves + f,1);

% Keep track of states
nStates_G = size(G_d.A,1);
x_G = zeros(nStates_G,kFinalWind+1);

% Set initial condition
x0_G = zeros(nStates_G,1);
x_G(:,1) = x0_G;

% Keep track of input and output sequences
uSeqFB = zeros(1,kFinalWind);
uSeqFF = zeros(1,kFinalWind);
uTotal = zeros(1,kFinalWind);
out = zeros(nOutputs,kFinalWind);

%%
% % Wind:
% % Turbulent wind
load('inputData\turbWind_16mps.mat') 
% v = windData;
v = [windData; windData; windData];
v = v-16; % center around linearization point
v = v./v0hat_max;

% EOG
% load('inputData\eog_16mps.mat','Wind1VelX','Time')
% v = interp1(Time,Wind1VelX,0:Ts_wind:Ts_wind*(600-1))'; % resample with current sampling period
% v = v-16;
% v = [zeros(500,1); v ; zeros(2000,1)];
% v = v./v0hat_max;
Mp = M_pitch(length(u_Mp)+1:end)./MpitchHat_max;
MpDeePC = Mp(1:20:end);
% v = zeros(20000,1);
% Mp = zeros(5000,1);


%% Solve the constrained optimization problem
ivFlag = 1;
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

ffFlag = 1; % enables/disables wave DeePC controller
fbFlag = 1; % enables/disables PI controller

for idxWind=1:kFinalWind
    disp('Iteration: ')
    disp(idxWind)

    if ffFlag == 1
        if mod(idxWind-1,Ts_ratio) == 0
            % Reference trajectory: zero reference, disturbance rejection
            rfWaves = refWaves((idxWaves-1)*nOutputs+1:(idxWaves-1)*nOutputs+nOutputs*f);

            % Wave preview
            data.wf = MpDeePC(idxWaves:idxWaves + controlParams.f - 1);
            uStarFF = deepc2(data,rfWaves,controlParams,method,ivFlag,1);
            uSeqFF(:,idxWind) = uStarFF(1);

            idxWaves = idxWaves + 1;
        else
            uSeqFF(:,idxWind) = uSeqFF(:,idxWind-1);
        end
    else
        uSeqFF(:,idxWind) = 0;
    end

    if ffFlag == 1 && fbFlag == 1
        if idxWind>1
            measOutput = out(:,idxWind-1);
            prevPitCom = uSeqFB(:,idxWind-1);
            prevIntErr = intErr;
        end
            %%%%%%%%%%% measOutput averaging %%%%%%%%%%%
            % avgWindowPoints = 1;
            % if idxWind < avgWindowPoints
            %     avgWindow = [zeros(avgWindowPoints-idxWind,1); out(:,1:idxWind)'];
            % else
            %     avgWindow = out(:,idxWind-avgWindowPoints+1:idxWind);
            % end
            % avgMeasOutput = mean(avgWindow);
            % 
            % [uSeqFB(:,idxWind),intErr] = getPIcom(refWind(idxWind),avgMeasOutput,prevPitCom,prevIntErr,Ts_wind);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           %%%%%%%%%%%%%%%%%%% LPF %%%%%%%%%%%%%%%%%%%%
            if idxWind > 25
                filteredOutput = lowpass(out(:,idxWind-25:idxWind),3e-3,Ts_wind);
            else
                filteredOutput = lowpass(out(:,1:idxWind),3e-3,Ts_wind);
            end

            [uSeqFB(:,idxWind),intErr] = getPIcom(refWind(idxWind),filteredOutput(end),prevPitCom,prevIntErr,Ts_wind);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           % [uSeqFB(:,idxWind),intErr] = getPIcom(refWind(idxWind),measOutput,prevPitCom,prevIntErr,Ts_wind);
        
    elseif ffFlag == 0 && fbFlag == 1
        if idxWind>1
            measOutput = out(:,idxWind-1);
            prevPitCom = uSeqFB(:,idxWind-1);
            prevIntErr = intErr;
        end
            
           %%%%%%%%%%%%%%%%%%% LPF %%%%%%%%%%%%%%%%%%%%
            if idxWind > 25
                filteredOutput = lowpass(out(:,idxWind-25:idxWind),3e-3,Ts_wind);
            else
                filteredOutput = lowpass(out(:,1:idxWind),3e-3,Ts_wind);
            end

            [uSeqFB(:,idxWind),intErr] = getPIcom(refWind(idxWind),filteredOutput(end),prevPitCom,prevIntErr,Ts_wind);
    end

    uTotal(:,idxWind) = uSeqFB(:,idxWind) + uSeqFF(:,idxWind);
    uTotal(:,idxWind) = min( max(uTotal(:,idxWind), deg2rad(-10)/uhat_max ),  deg2rad(10)/uhat_max );

    u = [uTotal(:,idxWind);
        v(idxWind)
        Mp(idxWind)];

    % Apply optimal input, simulate output
    x_G(:,idxWind+1) = G_d.A*x_G(:,idxWind) + G_d.B*u;
    out(:,idxWind) = G_d.C*x_G(:,idxWind) + G_d.D*u + Std.*randn(size(out(:,idxWind)));    

    if ffFlag == 1
        if mod(idxWind-1,Ts_ratio) == 0
            data.uini = [data.uini(nInputs+1:end); uSeqFF(:,idxWind); uSeqFB(:,idxWind)];
            data.yini = [data.yini(nOutputs+1:end); out(:,idxWind)];
            if previewFlagWaves == 1
                data.wini = [data.wini(nDist+1:end); Mp(idxWind)];
            end
        end
    end
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
yline(rad2deg(controlParams.lbu*uhat_max),'r--','LineWidth',1)
yline(rad2deg(controlParams.ubu*uhat_max),'r--','LineWidth',1)
% xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input FB')
grid on
set(gcf,'Color','White')

% Control input sequence FF
figure
plot(tsimWind,rad2deg(uSeqFF.*uhat_max))
hold on
plot(tsimWind,rad2deg(uSeqFF.*uhat_max))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-15 15])
yline(rad2deg(controlParams.lbu*uhat_max),'r--','LineWidth',1)
yline(rad2deg(controlParams.ubu*uhat_max),'r--','LineWidth',1)
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
yline(rad2deg(controlParams.lbu*uhat_max),'r--','LineWidth',1)
yline(rad2deg(controlParams.ubu*uhat_max),'r--','LineWidth',1)
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