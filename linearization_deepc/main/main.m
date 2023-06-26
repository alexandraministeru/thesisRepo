% This script obtains a reduced state-space (inputs: horizontal wind speed
% and collective blade pitch angle, outputs: generator speed) from a WT
% linearization, gets I/O open-loop data, and closes the control loop with
% a DeePC controller. Only the blade pitch angle is considered a
% controllled input, while the horizontal wind speed is a disturbance.

%% Clean environment
clearvars;clc;close all;

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
addpath(genpath('D:\Program Files\MATLAB\R2023a\casadi'));
import casadi.*

%% Plot nonlinear model outputs
% outFile = '..\forLin\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
% plotChannels = {'RotSpeed','GenSpeed','BldPitch1'};
% PlotFASToutput(outFile,[],[],plotChannels,1)

%% Load linearization
% My linearization
% FilePath = "..\5MW_OC3Spar_DLL_WTurb_WavesIrr";
% MtlbTlbxPath = "D:\Master\TUD\Y2\Thesis\matlab\repo\linearization_deepc\fromAmr\matlab-toolbox-main";
% [LTIsys, MBC, matData, FAST_linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath);
% cd ..\main

% LTIsys from Amr
FilePath = "..\fromAmr\5MW_OC3Spar_DLL_WTurb_WavesIrr_16mps";
MtlbTlbxPath = "D:\Master\TUD\Y2\Thesis\matlab\repo\linearization_deepc\fromAmr\matlab-toolbox-main";
[LTIsys, MBC, matData, FAST_linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath);
cd ..\..\main

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

%% PI controller
% Kp = -5e-3;
% Ki = -1e-4;
% 
% C = pid(Kp,Ki);
% C.u = 'e'; C.y = 'u';
% G_theta_wg.u = 'u'; G_theta_wg.y = 'y1';
% G_theta_v.u = 'v'; G_theta_v.y = 'y2';
% 
% Sum1 = sumblk('e = r - y');
% Sum2 = sumblk('y = y1 + y2');
% CL = connect(G_theta_wg,G_theta_v,C,Sum1,Sum2,{'r','v'},{'y'});
% 
% t = 0:1:500;
% r = zeros(size(t));
% v = zeros(size(t));
% v(100:end) = 1;
% 
% u = [r;v];
% 
% x0 = zeros(size(CL.A,1),1);
% x = x0;
% y = zeros(size(t));
% for k = 1:length(t)
%     xdot = CL.A*x + CL.B*u(:,k);
%     y(k) = CL.C*x +  CL.D*u(:,k);
%     x = xdot;
% end

%% Discretize and get OL data
% Discretize system
Ts = 0.05;

% sys = ss(G_theta_wg);
% sys_d = c2d(sys,Ts);

sys_d = c2d(LTIsys_reduced,Ts);

% Time vector
t = 0:Ts:50;

% PRBS input for persistency of excitation
u_bladePitch = idinput(length(t),'PRBS',[0 1/10],[-0.2 0.2]); % in degrees
u_bladePitch = u_bladePitch*(pi/180); % in radians
u_windSpeed = zeros(size(u_bladePitch));

u = [u_windSpeed u_bladePitch];

figure
plot(t,u_bladePitch)
grid on
xlabel('Time (in s)')
ylabel('Blade pitch angle (in rad)')
title('PRBS input')
set(gcf,'Color','White')

% Simulate system and get open-loop I/O data set
y = lsim(sys_d,u,t); % zero initial condition

% Add noise on output
noiseAns = input('Add noise? y/n[n]: ','s');
if isempty(noiseAns) || strcmp(noiseAns,'n')
    noiseFlag = 0;
elseif strcmp(noiseAns,'y')
    noiseFlag = 1;
end

variance = 0.1;
variance = variance*noiseFlag;
y = y + variance.*randn(size(y));

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

%% Construct data matrices
% Past data
data.Up = constructHankelMat(u_bladePitch,i,p,Nbar);
data.Yp = constructHankelMat(y,i,p,Nbar);

% Future data
data.Uf = constructHankelMat(u_bladePitch,i+p,f,Nbar);
data.Yf = constructHankelMat(y,i+p,f,Nbar);

%% Set up control loop
kFinal = 200; % simulation steps
tsim = 0:Ts:Ts*(kFinal-1);

% Generate reference trajectory
ref = zeros(kFinal+f,1);
ref(100:end) = 100;

% Keep track of states
nStates = size(sys_d.A,1);
x = zeros(nStates,kFinal+1);

% Set initial condition
x0 = zeros(nStates,1); % zero initial condition
x(:,1) = x0;

% Keep track of input and output sequences
uSeq = zeros(kFinal,1);
out = zeros(kFinal,1); % HARDCODED OUTPUT SIZE

%% Simulate/add wind disturbance
v = zeros(kFinal,1);
% % Step in incoming wind speed
% v(100:end) = 1;

% Turbulent wind obtained from previous simulation
% load('turbWind_16mps.mat')
% v = windData(1:kFinal);

%% Solve the constrained optimization problem
% Past data for prediction
uini = constructHankelMat(u_bladePitch,i+N-p,p,1);
yini = constructHankelMat(y,i+N-p,p,1);

% Define weights
Q = 1e2*eye(f);
R = eye(f);

% Choose optimization method
method = input(['Optimization method: 1-fmincon+SQP, 2-quadprog, ' ...
    '3-casadi+nlp: ']);

% Use instrumental variables
ivAns = input('Use instrumental variables? y/n[n]: ','s');
if isempty(ivAns) || strcmp(ivAns,'n')
    ivFlag = 0;
elseif strcmp(ivAns,'y')
    ivFlag = 1;
end

% Control loop
for k=1:kFinal
    disp('Iteration: ')
    disp(k)

    % Reference trajectory
    rf = ref(k:k+f-1);

    % DeePC optimal control input
    uStar = deepc(data,uini,yini,N,p,f,rf,Q,R,method,ivFlag);
    uSeq(k) = uStar;

    u = [v(k);
        uStar];
    
    % Apply optimal input, simulate output
    x(:,k+1) = sys_d.A*x(:,k) + sys_d.B*u;
    out(k) = sys_d.C*x(:,k) + sys_d.D*u + variance.*randn(1);
    
    % Update past data with most recent I/O data
    uini = [uini(2:end); uStar];
    yini = [yini(2:end); out(k)];
end

%% Plotting
figure
plot(tsim,out)
xlabel('Time (in s)')
ylabel('Generator speed (in rpm)')
title('DeePC for reference tracking')
grid on
hold on
plot(tsim,ref(1:kFinal)) % reference
legend('Controlled output','Reference','Location','SouthEast')
set(gcf,'Color','White')

% Control input sequence
figure
plot(tsim,uSeq*(180/pi))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
title('Control input')
grid on
set(gcf,'Color','White')
