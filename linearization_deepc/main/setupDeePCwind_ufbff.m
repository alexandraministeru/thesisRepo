function [dataWind,controlParamsWind] = setupDeePCwind(G,Ts)
% Scaling factors
uhat_max = 10*(pi/180); % Maximum expected input (rad)
v0hat_max = 0.1*16; % Maximum expected wind disturbance (m/s)
% Discretize system
G_d = c2d(G,Ts);

% Time vector
simTime = Ts*1000;
t = 0:Ts:simTime;

% Generate PRBS input for persistency of excitation
u_hat_bladePitch = idinput(length(t),'PRBS',[0 1/10],[-2 2]); % in degrees
u_hat_bladePitch = u_hat_bladePitch*(pi/180); % in radians

% Wind disturbance input
windData = 0.5.*randn(3000,1);
u_hat_windSpeed = windData(1:size(u_hat_bladePitch,1));

% No Waves
u_hat_Mp = zeros(length(t),1);

% Scale inputs
u_bladePitch = u_hat_bladePitch./uhat_max;
u_windSpeed = u_hat_windSpeed./v0hat_max;
u_Mp = u_hat_Mp;

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

% Generate PRBS input for persistency of excitation
u_hat_bladePitch_FF = idinput(length(t),'PRBS',[0 1/10],[-0.1 0.1]); % in degrees
u_hat_bladePitch_FF = u_hat_bladePitch_FF*(pi/180); % in radians
u_bladePitch_FF = u_hat_bladePitch_FF./uhat_max;


dataWind.u_bladePitch_FB = u_bladePitch - u_bladePitch_FF;
dataWind.u_bladePitch_FF = u_bladePitch_FF;
dataWind.u_bladePitch = u_bladePitch;
dataWind.u_windSpeed = u_windSpeed;
dataWind.y = y;

%% DeePC parameters
N = 500; % lenght of data set
p = 40; % past data window
f = 20; % prediction window
Nbar = N-p-f+1;
i = 1;

controlParamsWind.N = N;
controlParamsWind.p = p;
controlParamsWind.f = f;

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
dataWind.Up = constructHankelMat(dataWind.u_bladePitch,i,p,Nbar);
dataWind.Yp = constructHankelMat(dataWind.y,i,p,Nbar);

% Future data
dataWind.Uf = constructHankelMat(dataWind.u_bladePitch,i+p,f,Nbar);
dataWind.Yf = constructHankelMat(dataWind.y,i+p,f,Nbar);

disturbMat = [u_bladePitch_FF u_windSpeed];

if previewFlag == 1    
    dataWind.Wp = constructHankelMat(disturbMat,i,p,Nbar); % past data
    dataWind.Wf = constructHankelMat(disturbMat,i+p,f,Nbar); % future data
else
    dataWind.Wp = [];
    dataWind.Wf = [];
end

dataWind.uini = constructHankelMat(dataWind.u_bladePitch,i+N-p,p,1);
dataWind.yini = constructHankelMat(dataWind.y,i+N-p,p,1);

if previewFlag == 1
    dataWind.wini = constructHankelMat(disturbMat,i+N-p,p,1);
else
    dataWind.wini = [];
end

%% Set control weights
% Define control weights:
% weightOutputs diagonal matrix of size l-by-l, where l is the number of
% output channels and the n-th element on the diagonal represents the
% weight for the corresponding n-th output

weightOutputs = 5e2*diag(1); %5e2
controlParamsWind.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
weightInputs= 1*diag(1);
controlParamsWind.R = kron(eye(f),weightInputs);

% Choose input bounds
controlParamsWind.lbu = -10*(pi/180);
controlParamsWind.lbu = controlParamsWind.lbu/uhat_max;
controlParamsWind.ubu = 10*(pi/180);
controlParamsWind.ubu = controlParamsWind.ubu/uhat_max;

% Input rate constraint
duDeg = 8; % deg/s
duRad = duDeg*(pi/180); % rad/s
duRad = duRad/uhat_max;
controlParamsWind.duf = duRad*Ts;

end



