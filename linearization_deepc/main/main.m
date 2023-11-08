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
addpath(genpath('functions'))
% Insert full path to casadi lib
addpath(genpath('D:\Program Files\MATLAB\R2023a\casadi'));
import casadi.*
% YALMIP and solvers
addpath(genpath('D:\Program Files\mosek\10.1\toolbox\r2017aom')) %MOSEK
addpath(genpath('D:\Program Files\MATLAB\R2023a\sedumi')) % sedumi
addpath(genpath('D:\Program Files\MATLAB\R2023a\sdpt3')) % sdpt3
addpath(genpath('D:\Program Files\MATLAB\R2023a\yalmip')) % yalmip

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
% save('inputData\linData.mat','VTK','FAST_linData','LTIsys','matData','MBC');

% load('inputData\linData.mat');
load('inputData\linDataWave.mat');
load('inputData\waveForces.mat');

%% Find index of blade pitch, gen speed and wind speed
inputChannelsList = MBC.DescCntrlInpt;
outputChannelsList = MBC.DescOutput;

inputChannels = {'ED Extended input: collective blade-pitch command, rad', ...
    'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s',...
    'ED Platform Y moment, node 1, Nm'};

% inputChannels = {'ED Extended input: collective blade-pitch command, rad', ...
%     'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s'};
% outputChannels = {'ED GenSpeed, (rpm)'};
outputChannels = {'ED RotSpeed, (rpm)'};

LTIsys_reduced = getReducedSS(MBC,LTIsys,inputChannels,outputChannels);

% % Set input, output and state names
% LTIsys_reduced.InputName = {'Collective blade pitch (\theta_c)', 'Horizontal wind speed (v_0)'};
% LTIsys_reduced.InputUnit = {'rad','m/s'};
% % LTIsys_reduced.OutputName = {'Generator speed (\Omega_g)'};
% LTIsys_reduced.OutputName = {'Rotor speed (\Omega)'};
% LTIsys_reduced.OutputUnit = {'rpm'};

%% Discretize and get OL data
% Discretize system
Ts = 0.05; % sampling period (in s)
% Ts = 1;
sys_d = c2d(LTIsys_reduced,Ts);

% Time vector
simTime = Ts*1000;
t = 0:Ts:simTime;

% % Generate PRBS input for persistency of excitation
u_bladePitch = idinput(length(t),'PRBS',[0 1/10],[-2 2]); % in degrees
u_bladePitch = u_bladePitch*(pi/180); % in radians
 
% % Generate horizontal wind speed disturbance input
windData = 0.5.*randn(3000,1);
u_windSpeed = windData(1:size(u_bladePitch,1));
 
% % Save input to file
% save('inputData\peinput.mat','u_bladePitch','u_windSpeed','windData');

% Load input from file
% load('inputData\peinput.mat')

% % Can replace wind disturbance:
% u_windSpeed = zeros(size(u_bladePitch)); % no wind disturbance

u = [u_bladePitch u_windSpeed zeros(size(u_windSpeed))];

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

% std = 0.1; % measurement noise standard deviation
% std = 0.01; % measurement noise standard deviation for rotor speed control
std = 5e-3;
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
p = 40; % past data window
f = 20; % prediction window
Nbar = N-p-f+1;
i = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ts = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p = 20; % past data window
% f = 5; % prediction window
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
kFinal = 600; % simulation steps
% kFinal = 600*20;
tsim = 0:Ts:Ts*(kFinal-1);
nInputs = size(data.Up,1)/p;
nOutputs = size(data.Yp,1)/p;
if previewFlag == 1
    nDist = size(data.Wp,1)/p;
end

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
nStates = size(sys_d.A,1);
x = zeros(nStates,kFinal+1);

% Set initial condition
x0 = zeros(nStates,1); % zero initial condition
x(:,1) = x0;

% Keep track of input and output sequences
uSeq = zeros(nInputs,kFinal);
out = zeros(nOutputs,kFinal);

%% Wind disturbance for CL
% No wind disturbance
% v = zeros(kFinal+f,1);

% White noise
% v = windData(size(u_windSpeed,1)+1:end); 

% % Extreme operating gust
% load('inputData\eog_16mps.mat','Wind1VelX','Time')
% v = interp1(Time,Wind1VelX,tsim)'; % resample with current sampling period
% v = v-16; % center around linearization point
% v = [v;zeros(f,1)];

% Turbulent wind
load('inputData\turbWind_16mps.mat') %turbulent wind obtained from a previous FAST simulation
v = windData;
v = v-16; % center around linearization point

Mp = M_pitch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ts = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % load('inputData\eog_16mps.mat','Wind1VelX','Time')
% % v = Wind1VelX(1:125:end);
% % v = v-16;
% % v = [zeros(100,1); v; zeros(kFinal,1)];
% 
% % load('inputData\turbWind_16mps_long.mat')
% % % v = turbWind(1:20:end);
% % v = turbWind;
% % v = v-16; % center around linearization point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve the constrained optimization problem

% % Use instrumental variables
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
% weightOutputs diagonal matrix of size l-by-l, where l is the number of 
% output channels and the n-th element on the diagonal represents the 
% weight for the corresponding n-th output

if ivFlag == 1  
    % weightOutputs = 5e-3*diag(1); % genspeed
    weightOutputs = 5e-1*diag(1); % rotspeed
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

% Input rate constraint
duDeg = 8; % deg/s
duRad = 8*(pi/180); % rad/s
controlParams.duf = duRad*Ts;

% Choose optimization method
% method = input(['Optimization method: 1 - QP, ' '2 - SDP, ' '3 - NLP: ']);
method = 1;

%% Control loop
% Past data for prediction
data.uini = constructHankelMat(u_bladePitch,i+N-p,p,1);
data.yini = constructHankelMat(y,i+N-p,p,1);
if previewFlag == 1
    data.wini = constructHankelMat(u_windSpeed,i+N-p,p,1);
else
    data.wini = [];
end

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

    u = [uSeq(:,k);

        v(k); Mp(k)];
    
    % Apply optimal input, simulate output
    x(:,k+1) = sys_d.A*x(:,k) + sys_d.B*u;
    out(:,k) = sys_d.C*x(:,k) + sys_d.D*u + std.*randn(size(out(:,k)));
    
   % Update past data with most recent I/O data
    data.uini = [data.uini(nInputs+1:end); uStar];
    data.yini = [data.yini(nOutputs+1:end); out(:,k)];
    if previewFlag == 1
        data.wini = [data.wini(nDist+1:end); v(k)];
    end
end

%% Plotting
stepIdxs = find(ref);

% Controlled output
figure
plot(tsim,out)
xlabel('Time (in s)')   
ylabel([LTIsys_reduced.OutputName ' (in ' LTIsys_reduced.OutputUnit ')'])
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
plot(tsim,rad2deg(uSeq))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-15 15])
yline(rad2deg(controlParams.lbu),'r--','LineWidth',1)
yline(rad2deg(controlParams.ubu),'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input')
grid on
set(gcf,'Color','White')

% Control input rate
figure
plot(tsim(1:end-1),diff(rad2deg(uSeq))*(1/Ts))
xlabel('Time (in s)')
ylabel('\theta_c rate (in deg/s)')
yline(duDeg,'r--','LineWidth',1)
yline(-duDeg,'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
ylim([-15 15])
title('Control input rate')
grid on
set(gcf,'Color','White')

% Wind disturbance
figure
plot(tsim,v(1:length(tsim)))
xlabel('Time (in s)')
ylabel('Wind speed (in m/s)')
title('Horizontal wind speed around 16mps linearization point')
grid on
set(gcf,'Color','White')

%% Save data
% |   simType   |  description  | optMethod | scaled |     tuning     |
% |-------------|---------------|-----------|--------|----------------|
% |    wind     | noPreview     |    QP     |   []   |       []       |
% |    wave     |(wind)preview  |    SDP    | scaled |tuning*tunedVar*|
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
    simType = 'turb';  
    % simType = 'turbWind';
    inName = '\theta_c (in deg)';
    % outName = 'Generator speed (in rpm)';
    outName = 'Rotor speed (in rpm)';

    scaledFlag = 0;
    tuning = 0;
    tunedVar = 'Q';
    ext = '.mat';

    % Specify type of simulation
    fileName = simType;

    % Set description
    if previewFlag == 1
        if exist('disturbMat','var') == 0
            description = 'preview';
        elseif size(disturbMat,2) == 2
            description = 'Fsg_Mp_preview';
        elseif isequal(distrbMat,u_fsg)
             description = 'Fsg_preview';
        else
             description = 'Mp_preview';
        end  
        fileName = [fileName '_' description];
    else
        description = 'noPreview';
        fileName = [fileName '_' description];
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
    while isfile(['masterMeeting\' fileName ext])
        if isstrprop(fileName(end),'digit')
            currIdx = str2double(fileName(end));
            fileName = [fileName(1:end-1) num2str(currIdx+1)];
        else
            fileName = [fileName num2str(1)]; %#ok<AGROW>
        end
    end

    % Save variables
    if scaledFlag == 1
        save(['masterMeeting\' fileName ext],'inName','outName','tsim','kFinal','ref','Ts','controlParams', 'out', ...
            'uSeq','description','scaledFlag','uhat_max','ehat_max');
    else
        save(['masterMeeting\' fileName ext],'inName','outName','tsim','kFinal','ref','Ts','controlParams', 'out', ...
            'uSeq','description','scaledFlag');
    end
end


