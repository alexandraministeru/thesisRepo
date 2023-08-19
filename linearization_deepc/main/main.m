% This script obtains a reduced state-space (inputs: horizontal wind speed
% and collective blade pitch angle, outputs: generator speed) from a WT
% linearization, gets I/O open-loop data, and closes the control loop with
% a DeePC controller. Only the blade pitch angle is considered a
% controllled input, while the horizontal wind speed is a disturbance. The
% following parameters can be set during run time:
%
% 1. noiseFlag - selects whether or not measurement noise is added to the
% output of the system
% 2. previewFlag - selects whether or not preview disturbance information
% is provided to DeePC
% 3. ivFlag - selects whether or not instrumental variables are used within
% DeePC

%% Clean environment
clearvars;clc;close all;
rng('default')

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
% Insert full path to casadi lib
addpath(genpath('D:\Program Files\MATLAB\R2023a\casadi')); 
import casadi.*

%% Load linearization
% My linearization
% FilePath = "..\5MW_OC3Spar_DLL_WTurb_WavesIrr";
% MtlbTlbxPath = "D:\Master\TUD\Y2\Thesis\matlab\repo\linearization_deepc\fromAmr\matlab-toolbox-main";
% [LTIsys, MBC, matData, FAST_linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath);
% cd ..\main

% % LTIsys from Amr
% FilePath = "..\fromAmr\5MW_OC3Spar_DLL_WTurb_WavesIrr_16mps";
% MtlbTlbxPath = "D:\Master\TUD\Y2\Thesis\matlab\repo\linearization_deepc\fromAmr\matlab-toolbox-main";
% [LTIsys, MBC, matData, FAST_linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath);
% cd ..\..\main
% 
% save('linData.mat','VTK','FAST_linData','LTIsys','matData','MBC');

load('inputData\linData.mat');

%% Find index of blade pitch, gen speed and wind speed
inputChannelsList = MBC.DescCntrlInpt;
outputChannelsList = MBC.DescOutput;

inputChannels = {'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s',...
    'ED Extended input: collective blade-pitch command, rad'};
outputChannels = {'ED GenSpeed, (rpm)'};

nIn = length(inputChannels);
nOut = length(outputChannels);
nStates = length(LTIsys.A);

B_reduced = zeros(nStates,nIn);
C_reduced = zeros(nOut,nStates);
D_temp = zeros(length(LTIsys.C),nIn);
D_reduced = zeros(nOut,nIn);

for idx = 1:nIn
    % Find index of input channel
    id = find(ismember(inputChannelsList,inputChannels{idx}));
    B_reduced(:,idx) = MBC.AvgB(:,id);
    D_temp(:,idx) = MBC.AvgD(:,id);
end

for idx = 1:nOut
    % Find index of output channel
    id = find(ismember(outputChannelsList,outputChannels{idx}));
    C_reduced(idx,:) = MBC.AvgC(id,:);
    D_reduced(idx,:) = D_temp(id,:);
end

LTIsys_reduced = ss(MBC.AvgA,B_reduced,C_reduced,D_reduced);

%% Extract transfer functions
% G_theta_wg = tf(LTIsys(11, 9)); % blade pitch (rad) to gen speed (rpm)
% G_theta_v = tf(LTIsys(11, 1)); % wind speed (m/s) to gen speed (rpm)

%% Discretize and get OL data
% Discretize system
Ts = 0.05; % sampling period (in s)

% sys = ss(G_theta_wg);
% sys_d = c2d(sys,Ts);

sys_d = c2d(LTIsys_reduced,Ts);

% Time vector
t = 0:Ts:50;

% % Generate PRBS input for persistency of excitation
u_bladePitch = idinput(length(t),'PRBS',[0 1/10],[-2 2]); % in degrees
u_bladePitch = u_bladePitch*(pi/180); % in radians
% 
% % Generate horizontal wind speed disturbance input
windData = 0.5.*randn(3000,1);
u_windSpeed = windData(1:size(u_bladePitch,1));
% 
% % Save input to file
% save('inputData\peinput.mat','u_bladePitch','u_windSpeed','windData');

% Load input from file
% load('inputData\peinput.mat')

% % Can replace wind disturbance:
% u_windSpeed = zeros(size(u_bladePitch)); % no wind disturbance
% % or
% load('inputData\turbWind_16mps.mat') %turbulent wind obtained from a previous FAST simulation
% u_windSpeed = windData(1:size(u_bladePitch,1));
% u_windSpeed = u_windSpeed-16; % around linearization point

u = [u_windSpeed u_bladePitch];

figure
plot(t,u_bladePitch)
grid on
xlabel('Time (in s)')
ylabel('Blade pitch angle (in rad)')
title('PRBS input')
set(gcf,'Color','White')

figure
plot(t,u_windSpeed)
grid on
xlabel('Time (in s)')
ylabel('Horizontal wind speed (in m/s)')
title('Disturbance')
set(gcf,'Color','White')

% Simulate system and get open-loop I/O data set
y = lsim(sys_d,u,t); % zero initial condition, around linearization point

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

std = 0.1; % measurement noise standard deviation
std = std*noiseFlag;
y = y + std.*randn(size(y));

figure
plot(t,y);
grid on
xlabel('Time (in s)')
ylabel('Generator speed (in rpm)')
title('OL response to PE input - \omega_g around linearization point')
set(gcf,'Color','White')

% % Check if data is persistently exciting
% data = iddata(y,u,Ts);
% pexcit(data)

%% DeePC parameters
N = 500; % lenght of data set
p = 50; % past data window
f = 50; % prediction window
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

% Past data
data.Up = constructHankelMat(u_bladePitch,i,p,Nbar);
data.Yp = constructHankelMat(y,i,p,Nbar);

% Future data
data.Uf = constructHankelMat(u_bladePitch,i+p,f,Nbar);
data.Yf = constructHankelMat(y,i+p,f,Nbar);

if previewFlag == 1
    data.Wp = constructHankelMat(u_windSpeed,i,p,Nbar); % past data
    data.Wf = constructHankelMat(u_windSpeed,i+p,f,Nbar); % future data
else
    data.Wp = [];
    data.Wf = [];
end

%% Set up control loop
kFinal = 300; % simulation steps
tsim = 0:Ts:Ts*(kFinal-1);
nInputs = size(data.Up,1)/p;
nOutputs = size(data.Yp,1)/p;

% Generate reference trajectory
ref = zeros(kFinal+f,1);
ref(200:end) = 100; % step in reference

% %MIMO:
% ref_y1 = zeros(kFinal+f,1)';
% ref_y2 = zeros(kFinal+f,1)';
% ref_ym = zeros(kFinal+f,1)';
% ref = [ref_y1; ref_y2; ref_yk];
% ref = reshape(ref,[],1);

% Keep track of states
nStates = size(sys_d.A,1);
x = zeros(nStates,kFinal+1);

% Set initial condition
x0 = zeros(nStates,1); % zero initial condition
x(:,1) = x0;

% Keep track of input and output sequences
uSeq = zeros(nInputs,kFinal);
out = zeros(nOutputs,kFinal);

%% Wind disturbance for CL
v = windData(size(u_windSpeed,1)+1:end); 
% v = zeros(kFinal+f,1);

%% Solve the constrained optimization problem
% Past data for prediction
data.uini = constructHankelMat(u_bladePitch,i+N-p,p,1);
data.yini = constructHankelMat(y,i+N-p,p,1);

if previewFlag == 1
    data.wini = constructHankelMat(u_windSpeed,i+N-p,p,1);
else
    data.wini = [];
end

% Use instrumental variables
ivAns = input('Use instrumental variables? y/n[y]: ','s');

while not(isempty(ivAns)) && not(strcmp(ivAns,'y')) && ...
        not(strcmp(ivAns,'n'))    
    disp('Invalid input.')
    ivAns = input('Use instrumental variables? y/n[y]: ','s');
end

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
    weightOutputs = 5e-3*diag(1);
else
    weightOutputs = 1e4*diag(1);
end
controlParams.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of 
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
weightInputs= diag(1); 
controlParams.R = kron(eye(f),weightInputs);

% Choose input bounds
controlParams.lbu = -10*(pi/180);
controlParams.ubu = 10*(pi/180);

% Choose optimization method
method = input(['Optimization method: 1-quadprog, ' '2-casadi+nlp: ']);

% Control loop
for k=1:kFinal
    disp('Iteration: ')
    disp(k)

    % Reference trajectory
    rf = ref((k-1)*nOutputs+1:(k-1)*nOutputs+nOutputs*f);

    % Wind preview
    if previewFlag == 1
        data.wf = v(k:k+f-1);
    else
        data.wf = [];
    end

    % DeePC optimal control input
    uStar = deepc(data,rf,controlParams,method,ivFlag,previewFlag);
    uSeq(:,k) = uStar;

    u = [v(k);
        uStar];
    
    % Apply optimal input, simulate output
    x(:,k+1) = sys_d.A*x(:,k) + sys_d.B*u;
    out(:,k) = sys_d.C*x(:,k) + sys_d.D*u + std.*randn(size(out(:,k)));
    
   % Update past data with most recent I/O data
    data.uini = [data.uini(nInputs+1:end); uStar];
    data.yini = [data.yini(nOutputs+1:end); out(:,k)];
    if previewFlag == 1
        data.wini = [data.wini(nOutputs+1:end); v(k)];
    end
end

%% Plotting
stepIdxs = find(ref);

% Controlled output
figure
plot(tsim,out)
xlabel('Time (in s)')
ylabel('Generator speed (in rpm)')
title('DeePC for reference tracking')
grid on
hold on
plot(tsim,ref(1:kFinal)) % reference
xline(Ts*f,'k--','Future window size')
xline(Ts*stepIdxs(1),'k--','Reference step')
legend('Controlled output','Reference','Location','SouthEast')
set(gcf,'Color','White')

% Control input sequence
figure
plot(tsim,uSeq*(180/pi))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-15 15])
yline(-10,'r--','LineWidth',1)
yline(10,'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input')
grid on
set(gcf,'Color','White')

% plotCompare
