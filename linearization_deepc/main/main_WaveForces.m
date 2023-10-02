% DeePC with wave forces preview.
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
outFile = 'D:\Master\TUD\Y2\Thesis\matlab\fromAmr\FF_OpenFAST\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
[data, channels, units, headers] = ReadFASTtext(outFile);
plotChannels = {'RotSpeed','GenSpeed','BldPitch1','GenTq','GenPwr','Wave1Elev','B1WvsFxi','B1WvsMyi'};
PlotFASToutput(outFile,[],[],plotChannels,1)

F_surge = data(:,ismember(channels,'B1WvsFxi'));
M_pitch = data(:,ismember(channels,'B1WvsMyi'));

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

inputChannels = {'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s',...     
    'ED Extended input: collective blade-pitch command, rad', ...
    'ED Platform X force, node 1, N', ...
    'ED Platform Y moment, node 1, Nm'};
    % %, ... 
    %'ED Platform Y moment, node 1, Nm'};  

outputChannels = {'ED GenSpeed, (rpm)'};
    %,'ED RotSpeed, (rpm)', ...
    % 'HD B1WvsFxi, (N)', ...
    % 'HD B1WvsMyi, (N-m)'};

LTIsys_reduced = getReducedSS(MBC,LTIsys,inputChannels,outputChannels);

% save('inputData\waveLTIsys_reduced.mat','LTIsys_reduced');
% load('inputData\waveLTIsys_reduced.mat')

%% Discretize and get OL data
% Discretize system
Ts = 0.05; % sampling period (in s)
sys_d = c2d(LTIsys_reduced,Ts);

% Time vector
t = 0:Ts:50;

% Generate PRBS input for persistency of excitation
u_bladePitch = idinput(length(t),'PRBS',[0 1/10],[-2 2]); % in degrees
u_bladePitch = u_bladePitch*(pi/180); % in radians

% % Wind input: setady wind
% u_windSpeed = zeros(size(u_bladePitch));

% % Generate horizontal wind speed disturbance input
windData = 0.5.*randn(3000,1);
u_windSpeed = windData(1:size(u_bladePitch,1));

% Wave disturbance input: surge force, pitch moment
% u_Fsg = F_surge(1:length(u_bladePitch));
% u_Mp = M_pitch(1:length(u_bladePitch));

u_Fsg = std(F_surge).*randn(length(t),1);
u_Mp = std(M_pitch).*randn(length(t),1);

% % Save input to file
% save('inputData\peinput.mat','u_bladePitch');

% Load input from file
% load('inputData\peinput.mat')

% u = [u_windSpeed u_bladePitch u_Fsg];
u = [u_windSpeed u_bladePitch u_Fsg u_Mp];

figure
plot(t,u_bladePitch)
grid on
xlabel('Time (in s)')
ylabel('Blade pitch angle (in rad)')
title('PRBS input')
set(gcf,'Color','White')
% 
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
plot(t,y(:,1));
grid on
xlabel('Time (in s)')
ylabel('Generator speed (in rpm)')
title('OL response to PE input - \omega_g around linearization point')
set(gcf,'Color','White')

% % Check if data is persistently exciting
data = iddata(y,u,Ts);
pexcit(data)

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

clear data

% Past data
data.Up = constructHankelMat(u_bladePitch,i,p,Nbar);
data.Yp = constructHankelMat(y,i,p,Nbar);

% Future data
data.Uf = constructHankelMat(u_bladePitch,i+p,f,Nbar);
data.Yf = constructHankelMat(y,i+p,f,Nbar);

disturbMat = [u_Fsg u_Mp];

if previewFlag == 1
    data.Wp = constructHankelMat(disturbMat,i,p,Nbar); % past data
    data.Wf = constructHankelMat(disturbMat,i+p,f,Nbar); % future data
else
    data.Wp = [];
    data.Wf = [];
end

%% Set up control loop
kFinal = 1000; % simulation steps
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
nStates = size(sys_d.A,1);
x = zeros(nStates,kFinal+1);

% Set initial condition
x0 = zeros(nStates,1); % zero initial condition
x(:,1) = x0;

% Keep track of input and output sequences
uSeq = zeros(nInputs,kFinal);
out = zeros(nOutputs,kFinal);

%% CL disturbances
% Turbulent wind
% load('inputData\turbWind_16mps.mat') %turbulent wind obtained from a previous FAST simulation
% v = windData;
% v = v-16; % center around linearization point

v = zeros(kFinal+f,1); % Steady wind

% % Extreme operating gust
% load('inputData\eog_16mps.mat','Wind1VelX','Time')
% v = interp1(Time,Wind1VelX,tsim)'; % resample with current sampling period
% v = v-16; % center around linearization point
% v = [v;zeros(f,1)];

% Fsg = F_surge(length(u_Fsg)+1:end); 
% Mp = M_pitch(length(u_Mp)+1:end);

Fsg = F_surge; 
Mp = M_pitch;

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
    weightOutputs = 1e-2*diag(1);
else
    weightOutputs = 1e2*diag(1);
end
controlParams.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of 
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
weightInputs= 1*diag(1); 
controlParams.R = kron(eye(f),weightInputs);

% Choose input bounds
controlParams.lbu = -10*(pi/180);
controlParams.ubu = 10*(pi/180);

% Input rate constraint
duDeg = 8; % deg/s
duRad = 8*(pi/180); % rad/s
controlParams.duf = duRad*Ts;

% Choose optimization method
method = input(['Optimization method: 1-quadprog, ' '2-casadi+nlp: ']);

%% Control loop
tic
for k=1:1000
    disp('Iteration: ')
    disp(k)

    % Reference trajectory
    rf = ref((k-1)*nOutputs+1:(k-1)*nOutputs+nOutputs*f);

    % Wind preview
    if previewFlag == 1        
        % data.wf = [Mp(k:k+f-1)];
        previewData = [Fsg(k:k+f-1) Mp(k:k+f-1)];
        data.wf = reshape(previewData',[],1);
    else
        data.wf = [];
    end

    % DeePC optimal control input
    uStar = deepc(data,rf,controlParams,method,ivFlag,previewFlag);
    uSeq(:,k) = uStar;

    u = [v(k);
        uStar;
        Fsg(k);
        Mp(k)];        
    
    % Apply optimal input, simulate output
    x(:,k+1) = sys_d.A*x(:,k) + sys_d.B*u;
    out(:,k) = sys_d.C*x(:,k) + sys_d.D*u + std.*randn(size(out(:,k)));
    
   % Update past data with most recent I/O data
    data.uini = [data.uini(nInputs+1:end); uStar];
    data.yini = [data.yini(nOutputs+1:end); out(:,k)];
    if previewFlag == 1
        data.wini = [data.wini(nDist+1:end); Fsg(k); Mp(k)];
    end
end
toc

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
% xline(Ts*stepIdxs(1),'k--','Reference step')
legend('Controlled output','Reference','Location','SouthEast')
set(gcf,'Color','White')

% Control input sequence
figure
plot(tsim,rad2deg(uSeq))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-15 15])
yline(controlParams.lbu,'r--','LineWidth',1)
yline(controlParams.ubu,'r--','LineWidth',1)
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

%% Save data
% |   simType   |  description  | optMethod | scaled |     tuning     |
% |-------------|---------------|-----------|--------|----------------|
% |    wind     | noPreview     |    QP     |   []   |       []       |
% |    wave     |(wind)preview  |    SDP    | scaled |tuning*tunedVar*|
% |             | Fsg_preview   |    NLP    |        |                |
% |             | Mp_preview    |           |        |                |
% |             | Fsg_Mp_preview|           |        |                |
% ---------------------------------------------------------------------

simType = 'wave';
inName = '\theta_c (in deg)';
outName = 'Generator speed (in rpm)';
ext = '.mat';

% Specify type of simulation
fileName = simType;

% Set description
if previewFlag == 1
    if size(disturbMat,2) == 2
        description = 'Fsg_Mp_preview';
    else
        if exist('u_fsg','var') == 0
            description = 'preview';
        else
            if isequal(distrbMat,u_fsg)
                description = 'Fsg_preview';
            else
                description = 'Mp_preview';
            end
        end
    end
    fileName = [fileName '_' description];
else
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
    save(['outputData\' fileName ext],'inName','outName','tsim','ref','Ts','controlParams', 'out', ...
        'uSeq','description','scaledFlag','uhat_max','ehat_max');
else
    save(['outputData\' fileName ext],'inName','outName','tsim','ref','Ts','controlParams', 'out', ...
        'uSeq','description','scaledFlag');
end
