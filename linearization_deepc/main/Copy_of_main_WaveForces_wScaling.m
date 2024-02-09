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
% % plotChannels = {'RotSpeed','GenSpeed','BldPitch1','GenTq','GenPwr','Wave1Elev','B1WvsFxi','B1WvsMyi'};
% % PlotFASToutput(outFileOP,[],[],plotChannels,1)
%
% F_surge = data(:,ismember(channels,'B1WvsFxi'));
% M_pitch = data(:,ismember(channels,'B1WvsMyi'));
% save('inputData\waveForces.mat','F_surge','M_pitch');

load('inputData\waveForces.mat');
M_pitch = M_pitch(1:20:end); % 1/0.05

%% Load linearization
% FilePath = "..\5MW_OC3Spar_DLL_WTurb_WavesIrr";
% MtlbTlbxPath = "D:\Master\TUD\Y2\Thesis\matlab\repo\linearization_deepc\fromAmr\matlab-toolbox-main";
% [LTIsys, MBC, matData, FAST_linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath);
% cd ..\main
% save('inputData\linDataWave.mat','VTK','FAST_linData','LTIsys','matData','MBC');

load('inputData\linDataWave.mat');

%% Find index of blade pitch, gen speed and wind speed
inputChannelsList = MBC.DescCntrlInpt;
outputChannelsList = MBC.DescOutput;

inputChannels = {'ED Extended input: collective blade-pitch command, rad', ...
    'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s',...
    'ED Platform Y moment, node 1, Nm'};
    % 'ED Platform X force, node 1, N', ...

% outputChannels = {'ED GenSpeed, (rpm)'};
outputChannels = {'ED RotSpeed, (rpm)'};

LTIsys_reduced = getReducedSS(MBC,LTIsys,inputChannels,outputChannels);

% Set input, output and state names
LTIsys_reduced.InputName = {'Collective blade pitch (\theta_c)', ...
    'Horizontal wind speed (v_0)', ...
    'Wave pitch moment (M_{pitch})'}; 
% 'Wave surge force (F_{surge})',

LTIsys_reduced.InputUnit = {'rad','m/s','Nm'};

% LTIsys_reduced.OutputName = {'Generator speed (\Omega_g)'};
LTIsys_reduced.OutputName = {'Rotor speed (\Omega)'};

LTIsys_reduced.OutputUnit = {'rpm'};
% 
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
% save('inputData\waveLTIsys_reduced_rotSpeed.mat','LTIsys_reduced');

load('inputData\waveLTIsys_reduced_rotSpeed.mat')
% load('inputData\waveLTIsys_reduced.mat')

G_hat = LTIsys_reduced;

%% Scale transfer functions(G_hat, Gd_hat - unscaled; G, Gd scaled)
% Scaling factors
uhat_max = 10*(pi/180); % Maximum expected input (rad)
v0hat_max = 0.1*16; % Maximum expected wind disturbance (m/s)
FsurgHat_max = max(F_surge); % Maximum expected wave surge force (N)
MpitchHat_max = max(M_pitch); % Maximum expected wave pitch moment (Nm)
% ehat_max = 0.1* 1173.69; % Maximum expected generator speed error (10% around linearization OP) (rpm)
ehat_max = 0.1* 12.1; % Maximum expected generator speed error (10% around linearization OP) (rpm)

% Du = diag([uhat_max v0hat_max FsurgHat_max MpitchHat_max]);
Du = diag([uhat_max v0hat_max MpitchHat_max]);
Dy = ehat_max;

% Scaled transfer functions
G = Dy\G_hat*Du;

%% Discretize and get OL data
% Discretize system
Ts = 0.05; % sampling period (in s)
% Ts = 1;
G_d = c2d(G,Ts);

% Time vector
simTime = Ts*1000;
t = 0:Ts:simTime;

% Generate PRBS input for persistency of excitation
u_hat_bladePitch = idinput(length(t),'PRBS',[0 1/10],[-2 2]); % in degrees
u_hat_bladePitch = u_hat_bladePitch*(pi/180); % in radians

% Wind input: setady wind
% u_hat_windSpeed = zeros(size(u_hat_bladePitch));
windData = 0.5.*randn(3000,1);
u_hat_windSpeed = windData(1:size(u_hat_bladePitch,1));

% Wave disturbance input: surge force, pitch moment
u_hat_Fsg = std(F_surge).*randn(length(t),1);
% u_hat_Mp = std(M_pitch).*randn(length(t),1);
u_hat_Mp = zeros(length(u_hat_windSpeed),1);

% figure
% plot(t,u_bladePitch)
% grid on
% xlabel('Time (in s)')
% ylabel('Blade pitch angle (in rad)')
% title('PRBS input')
% set(gcf,'Color','White')

% figure
% plot(t,u_Fsg)
% grid on
% xlabel('Time (in s)')
% ylabel('Wave forces in surge direction (in N)')
% title('WvsFx')
% set(gcf,'Color','White')

% figure
% plot(t,u_Mp)
% grid on
% xlabel('Time (in s)')
% ylabel('Wave pitch moment on platform (in N-m)')
% title('WvsMy')
% set(gcf,'Color','White')

% Scale inputs
u_bladePitch = u_hat_bladePitch./uhat_max;
u_windSpeed = u_hat_windSpeed./v0hat_max;
u_Fsg = u_hat_Fsg./FsurgHat_max;
u_Mp = u_hat_Mp./MpitchHat_max;

% Simulate system and get open-loop I/O data set
y = lsim(G_d,[u_bladePitch u_windSpeed u_Mp],t); % zero initial condition, around linearization point

% Add measurement noise on output
noiseAns = input('Add measurement noise? y/n[y]: ','s');

while not(isempty(noiseAns)) && not(strcmp(noiseAns,'y')) && ...
        not(strcmp(noiseAns,'n'))
    disp('Invalid input.')
    noiseAns = input('Add measurement noise? y/n[y]: ','s');
end

if isempty(noiseAns) || strcmp(noiseAns,'y')
    noiseFlag = 1;
elseif strcmp(noiseAns,'n')
    noiseFlag = 0;
end

% std = 0.01;
Std = 5e-5; % measurement noise standard deviation
Std = Std*noiseFlag;
y = y + Std.*randn(size(y));

figure
plot(t,y(:,1))
grid on
xlabel('Time (in s)')
ylabel('Rotor speed, scaled (-)')
% ylabel('Generator speed (in rpm)')
title('OL response to PE input - \omega_g around linearization point')
set(gcf,'Color','White')

%% DeePC parameters
N = 500; % lenght of data set
p = 40; % past data window
f = 20; % prediction window
Nbar = N-p-f+1;
i = 1;

controlParams.N = N;
controlParams.p = p;
controlParams.f = f;

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

clear data

% Past data
data.Up = constructHankelMat(u_bladePitch,i,p,Nbar);
data.Yp = constructHankelMat(y,i,p,Nbar);

% Future data
data.Uf = constructHankelMat(u_bladePitch,i+p,f,Nbar);
data.Yf = constructHankelMat(y,i+p,f,Nbar);

% disturbMat = [u_Fsg u_Mp];
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
kFinal = 400; % simulation steps
tsim = 0:Ts:Ts*(kFinal-1);
nInputs = size(data.Up,1)/p;
nOutputs = size(data.Yp,1)/p;
nDist = size(data.Wp,1)/p;

% Generate reference trajectory
ref = zeros(kFinal+f,1);
% ref(200:end) = 100; % step in reference

% %MIMO:
% ref_y1 = zeros(kFinal+f,1)';
% ref_y2 = zeros(kFinal+f,1)';
% ref_ym = zeros(kFinal+f,1)';
% ref = [ref_y1; ref_y2; ref_yk];
% ref = reshape(ref,[],1);

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
% TurbWind
load('inputData\turbWind_16mps_long.mat')
v = turbWind(1:20:end);
v = v-16; % center around linearization point
v = v./v0hat_max;

% % EOG
% load('inputData\eog_16mps.mat','Wind1VelX','Time')
% v = Wind1VelX(1:125:end);
% v = v-16;
% v = [zeros(100,1); v; zeros(kFinal,1)];
% v = v./v0hat_max;

% v = zeros(kFinal+f,1)/v0hat_max; % Steady wind
% Fsg = F_surge./FsurgHat_max;
% Mp = M_pitch./MpitchHat_max ;
Mp = zeros(1000,1);

%% Solve the constrained optimization problem
% Past data for prediction
data.uini = constructHankelMat(u_bladePitch,i+N-p,p,1);
data.yini = constructHankelMat(y,i+N-p,p,1);

if previewFlag == 1
    data.wini = constructHankelMat(disturbMat,i+N-p,p,1);
else
    data.wini = [];
end

% Use instrumental variables
% ivAns = input('Use instrumental variables? y/n[y]: ','s');
% 
% while not(isempty(ivAns)) && not(strcmp(ivAns,'y')) && ...
%         not(strcmp(ivAns,'n'))
%     disp('Invalid input.')
%     ivAns = input('Use instrumental variables? y/n[y]: ','s');
% end
ivAns = 'y';

if isempty(ivAns) || strcmp(ivAns,'y')
    ivFlag = 1;
elseif strcmp(ivAns,'n')
    ivFlag = 0;
end

% Define control weights:
%
% weightOutputs diagonal matrix of size l-by-l, where l is the number of
% output channels and the n-th element on the diagonal represents the
% weight for the corresponding n-th output

if ivFlag == 1
    weightOutputs = 0.5*diag(1);
else
    weightOutputs = 2e2*diag(1);
end

controlParams.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
weightInputs= 1*diag(1);
controlParams.R = kron(eye(f),weightInputs);

% Choose input bounds
controlParams.lbu = -10*(pi/180);
controlParams.lbu = controlParams.lbu/uhat_max;
controlParams.ubu = 10*(pi/180);
controlParams.ubu = controlParams.ubu/uhat_max;

% Input rate constraint
duDeg = 8; % deg/s
duRad = duDeg*(pi/180); % rad/s
duRad = duRad/uhat_max;
controlParams.duf = duRad*Ts;

% Choose optimization method
% method = input(['Optimization method: 1 - QP, ' '2 - SDP, ' '3 - NLP: ']);
method = 1;

%% Control loop
tic
for k=1:kFinal
    disp('Iteration: ')
    disp(k)

    % Reference trajectory
    rf = ref((k-1)*nOutputs+1:(k-1)*nOutputs+nOutputs*f);

    % Wind preview
    if previewFlag == 1
        data.wf = v(k:k+f-1);
        % previewData = [v(k:k+f-1) Mp(k:k+f-1)];
        % data.wf = reshape(previewData',[],1);
    else
        data.wf = [];
    end

    % DeePC optimal control input
    uStar = deepc(data,rf,controlParams,method,ivFlag,previewFlag);
    uSeq(:,k) = uStar;

    u = [uStar;
        v(k);
        Mp(k)];

    % Apply optimal input, simulate output
    x_G(:,k+1) = G_d.A*x_G(:,k) + G_d.B*u;
    out(:,k) = G_d.C*x_G(:,k) + G_d.D*u + Std.*randn(size(out(:,k)));

    % Update past data with most recent I/O data
    data.uini = [data.uini(nInputs+1:end); uStar];
    data.yini = [data.yini(nOutputs+1:end); out(:,k)];
    if previewFlag == 1
        data.wini = [data.wini(nDist+1:end); v(k)];
    end
end
toc

%% Plotting
stepIdxs = find(ref);

% Controlled output
figure
plot(tsim,out.*ehat_max)
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
% ylabel('Generator speed (in rpm)')
title('DeePC for reference tracking')
grid on
hold on
plot(tsim,ref(1:kFinal)) % reference
xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
legend('Controlled output','Reference','Location','SouthEast')
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
plot(tsim(1:end-1), rad2deg((diff(uSeq).*uhat_max).*(1/Ts)))
xlabel('Time (in s)')
ylabel('\theta_c rate (in deg/s)')
yline(duDeg,'r--','LineWidth',1)
yline(-duDeg,'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
ylim([-15 15])
title('Control input rate')
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

saveAns = input('Save data? y/n[y]: ','s');
if isempty(saveAns) || strcmp(saveAns,'y')
    saveFlag = 1;
elseif strcmp(saveAns,'n')
    saveFlag = 0;
end

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
    if previewFlag == 1
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
    while isfile(['outputData\' fileName ext])
        if isstrprop(fileName(end),'digit')
            currIdx = str2double(fileName(end));
            fileName = [fileName(1:end-1) num2str(currIdx+1)];
        else
            fileName = [fileName num2str(1)]; %#ok<AGROW>
        end
    end

    % Save variables
    if scaledFlag == 1
        save(['outputData\' fileName ext],'inName','outName','tsim','kFinal','ref','Ts','controlParams', 'out', ...
            'uSeq','description','scaledFlag','uhat_max','ehat_max');
    else
        save(['outputData\' fileName ext],'inName','outName','tsim','kFinal','ref','Ts','controlParams', 'out', ...
            'uSeq','description','scaledFlag');
    end
end