function [dataWaves,controlParamsWaves] = setupDeePCwave(G,Ts,M_pitch)
% Scaling factors
uhat_max = 10*(pi/180); % Maximum expected input (rad)
v0hat_max = 0.1*16; % Maximum expected wind disturbance (m/s)
MpitchHat_max = max(M_pitch); % Maximum expected wave pitch moment (Nm)

% Discretize system
G_d = c2d(G,Ts);

% Time vector
simTime = Ts*1000;
t = 0:Ts:simTime;

% Generate PRBS input for persistency of excitation
u_hat_bladePitch = idinput(length(t),'PRBS',[0 1/10],[-2 2]); % in degrees
u_hat_bladePitch = u_hat_bladePitch*(pi/180); % in radians

% Steady wind
u_hat_windSpeed = zeros(length(t),1);

% Wave disturbance input: pitch moment
u_hat_Mp = std(M_pitch).*randn(length(t),1);

% Scale inputs
u_bladePitch = u_hat_bladePitch./uhat_max;
u_windSpeed = u_hat_windSpeed./v0hat_max;
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
title('OL response to PE input - \omega_g around linearization point')
set(gcf,'Color','White')

dataWaves.u_bladePitch = u_bladePitch;
dataWaves.u_Mp = u_Mp;
dataWaves.y = y;

%% DeePC parameters
N = 600; % lenght of data set
p = 20; % past data window
f = 20; % prediction window
Nbar = N-p-f+1;
i = 1;

controlParamsWaves.N = N;
controlParamsWaves.p = p;
controlParamsWaves.f = f;

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
dataWaves.Up = constructHankelMat(dataWaves.u_bladePitch,i,p,Nbar);
dataWaves.Yp = constructHankelMat(dataWaves.y,i,p,Nbar);

% Future data
dataWaves.Uf = constructHankelMat(dataWaves.u_bladePitch,i+p,f,Nbar);
dataWaves.Yf = constructHankelMat(dataWaves.y,i+p,f,Nbar);

if previewFlag == 1
    dataWaves.Wp = constructHankelMat(dataWaves.u_Mp,i,p,Nbar); % past data
    dataWaves.Wf = constructHankelMat(dataWaves.u_Mp,i+p,f,Nbar); % future data
else
    dataWaves.Wp = [];
    dataWaves.Wf = [];
end

dataWaves.uini = constructHankelMat(dataWaves.u_bladePitch,i+N-p,p,1);
dataWaves.yini = constructHankelMat(dataWaves.y,i+N-p,p,1);

if previewFlag == 1
    dataWaves.wini = constructHankelMat(dataWaves.u_Mp,i+N-p,p,1);
else
    dataWaves.wini = [];
end

%% Set control weights
% Define control weights:
% weightOutputs diagonal matrix of size l-by-l, where l is the number of
% output channels and the n-th element on the diagonal represents the
% weight for the corresponding n-th output

weightOutputs = 1e2*diag(1);
controlParamsWaves.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
weightInputs= 1*diag(1);
controlParamsWaves.R = kron(eye(f),weightInputs);

% Choose input bounds
controlParamsWaves.lbu = -deg2rad(10);
controlParamsWaves.lbu = controlParamsWaves.lbu/uhat_max;
controlParamsWaves.ubu = deg2rad(10);
controlParamsWaves.ubu = controlParamsWaves.ubu/uhat_max;

% Input rate constraint
duDeg = 8; % deg/s
duRad = duDeg*(pi/180); % rad/s
duRad = duRad/uhat_max;
controlParamsWaves.duf = duRad*Ts;

end



