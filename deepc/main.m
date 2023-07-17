clearvars;clc;close all;
addpath(genpath('D:\Program Files\MATLAB\R2023a\casadi'));
import casadi.*
rng('default')
%% Get I/O open-loop data
% Sampling time
Ts = 0.1;

% % Generate data
% [u,y,sys_d] = getIOData(Ts);
% save('iodata.mat','u','y','sys_d');

% Load data from file
load('iodata.mat')

% Add noise on output
noiseAns = input('Add noise? y/n[n]: ','s');
if isempty(noiseAns) || strcmp(noiseAns,'n')
    noiseFlag = 0;
elseif strcmp(noiseAns,'y')
    noiseFlag = 1;
end

std = 0.01;
std = std*noiseFlag;
y = y + std.*randn(size(y));

%% DeePC parameters
N = 200; % lenght of data set
p = 20; % past data window
f = 20; % prediction window
Nbar = N-p-f+1;
i = 1;

%% Construct data matrices
% Past data
data.Up = constructHankelMat(u,i,p,Nbar);
data.Yp = constructHankelMat(y,i,p,Nbar);

% Future data
data.Uf = constructHankelMat(u,i+p,f,Nbar);
data.Yf = constructHankelMat(y,i+p,f,Nbar);

%% Solve the constrained optimization problem
% Define control weights:
%
% weightOutputs diagonal matrix of size l-by-l, where l is the number of 
% output channels and the n-th element on the diagonal represents the 
% weight for the corresponding n-th input
weightOutputs = 5e3*diag(1); 
controlParams.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of 
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
weightInputs= diag(1); 
controlParams.R = kron(eye(f),weightInputs);

% Input bounds: lbu = [lbu1 lbu2 ... lbuk]', where k is the number of
% input channels; ubu is contstructed similarly
controlParams.lbu = -12;
controlParams.ubu = 12;

% Past data for prediction
uini = constructHankelMat(u,i+N-p,p,1);
yini = constructHankelMat(y,i+N-p,p,1);

kFinal = 100; % simulation steps
tsim = 0:Ts:Ts*(kFinal-1); % time vector
nStates = size(sys_d.A,1);
nInputs = size(sys_d.B,2);
nOutputs = size(sys_d.C,1);

% Generate reference trajectory - modify for MIMO
ref = zeros(kFinal+f,1);
ref(70:end) = 1;

% %MIMO:
% ref_y1 = zeros(kFinal+f,1)';
% ref_y2 = zeros(kFinal+f,1)';
% ref_ym = zeros(kFinal+f,1)';
% ref = [ref_y1; ref_y2; ref_yk];
% ref = reshape(ref,[],1);

% Keep track of states
x = zeros(nStates,kFinal+1);
x0 = zeros(nStates,1); % initial condition
% x0 = [0.5 0.2]; % initial condition
x(:,1) = x0;

% Keep track of output
out = zeros(nOutputs,kFinal);

% Control input sequence
uSeq = zeros(nInputs,kFinal);

% Choose optimization method
method = input(['Optimization method: 1-fmincon+SQP, 2-quadprog, ' ...
    '3-casadi+nlp, 4-casadi+QP: ']);

% Use instrumental variables
ivAns = input('Use instrumental variables? y/n[n]: ','s');
if isempty(ivAns) || strcmp(ivAns,'n')
    ivFlag = 0;
elseif strcmp(ivAns,'y')
    ivFlag = 1;
end

tic
% Control loop
for k=1:kFinal
    disp('Iteration: ')
    disp(k)

    % Reference trajectory
    rf = ref((k-1)*nOutputs+1:(k-1)*nOutputs+nOutputs*f);

    % DeePC optimal control input
    uStar = deepc(data,uini,yini,N,p,f,rf,controlParams,method,ivFlag);
    uSeq(:,k) = uStar;

    % Apply optimal input, simulate output including additive white noise
    x(:,k+1) = sys_d.A*x(:,k) + sys_d.B.*uStar;
    out(:,k) = sys_d.C*x(:,k) + sys_d.D.*uStar + std*randn(size(out(:,k)));

    % Update past data with most recent I/O data
    uini = [uini(nInputs+1:end); uStar];
    yini = [yini(nOutputs+1:end); out(:,k)];
end
toc

%% Plotting
figure
plot(tsim,out)
xlabel('Time (in s)')
ylabel('Displacement (in m)')
title('DeePC')
grid on
hold on
% Reference
plot(tsim,ref(1:kFinal))
xlabel('Time (in s)')
ylabel('Displacement (in m)')
legend('Controlled output','Reference','Location','SouthEast')
xline(f*Ts,'r--','Future window size')
set(gcf,'Color','White')

% Control input
figure
plot(tsim,uSeq)
xlabel('Time (in s)')
ylabel('Control input (in N)')
ylim([-14 14])
title('Control input')
grid on
set(gcf,'Color','White')
yline(-12,'r--','LineWidth',1)
yline(12,'r--','LineWidth',1)
xline(f*Ts,'k--','Future window size')

