clearvars;clc;close all;
addpath(genpath('D:\Program Files\MATLAB\R2023a\casadi'));
import casadi.*
%% Get I/O open-loop data
Ts = 0.1;
[u,y,sys_d] = getIOData(Ts);

% Add noise on output
noiseAns = input('Add noise? y/n[n]: ','s');
if isempty(noiseAns) || strcmp(noiseAns,'n')
    noiseFlag = 0;
elseif strcmp(noiseAns,'y')
    noiseFlag = 1;
end

variance = 0.01;
variance = variance*noiseFlag;
y = y + variance.*randn(size(y));

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
% Define weights
Q = 1e4*eye(f);
R = eye(f);

% Past data for prediction
uini = constructHankelMat(u,i+N-p,p,1);
yini = constructHankelMat(y,i+N-p,p,1);

kFinal = 100; % simulation steps
tsim = 0:Ts:Ts*(kFinal-1); % time vector
nStates = size(sys_d.A,1);
x0 = zeros(nStates,1); % initial condition
% x0 = [0.5 0.2]; % initial condition

% Generate reference trajectory
ref = ones(kFinal+f,1);

% Keep track of states
x = zeros(nStates,kFinal+1);
x(:,1) = x0;

% Keep track of output
out = zeros(kFinal,1);

% Control input sequence
uSeq = zeros(kFinal,1);

% Choose optimization method
method = input('Optimization method: 1-fmincon+SQP, 2-quadprog, 3-casadi+nlp: ');

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
    uStar = deepcFunc(data,uini,yini,N,p,f,rf,Q,R,method,ivFlag);
    uSeq(k) = uStar;

    % Apply optimal input, simulate output including additive white noise
    x(:,k+1) = sys_d.A*x(:,k) + sys_d.B.*uStar;
    out(k) = sys_d.C*x(:,k) + sys_d.D.*uStar + variance*randn(size(out(k)));

    % Update past data with most recent I/O data
    uini = [uini(2:end); uStar];
    yini = [yini(2:end); out(k)];
end

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
set(gcf,'Color','White')

% Control input
figure
plot(tsim,uSeq)
xlabel('Time (in s)')
ylabel('Control input (in N)')
title('Control input')
grid on
set(gcf,'Color','White')

