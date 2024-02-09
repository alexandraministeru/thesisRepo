%% Clean environment
clearvars;close all; clc
rng('default')
cd 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\MATLAB_Files\waveFF\'
addpath(genpath('functions'))

%% Set paths
% YALMIP and solvers
addpath(genpath('D:\Program Files\mosek\10.1\toolbox\r2017aom')) %MOSEK
addpath(genpath('D:\Program Files\MATLAB\R2023a\sedumi')) % sedumi
addpath(genpath('D:\Program Files\MATLAB\R2023a\sdpt3')) % sdpt3
addpath(genpath('D:\Program Files\MATLAB\R2023a\yalmip')) % yalmip

% Insert full path to casadi lib
addpath(genpath('D:\Program Files\MATLAB\R2023a\casadi'));
import casadi.*

%% Load library
addpath('..\..\');

% Load library
loadlibrary('D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\QBladeCE_2.0.6.dll','QBladeDLLFunctions.h','alias','QBladeDLL')
if not(libisloaded('QBladeDLL'))
    fprintf('Library not loaded')
end

m = libfunctions('QBladeDLL') ;

% Store library functions
if isempty(m)
    fprintf('Error')
end

% % Help function
% QBladeFunctionHelp('advanceController_at_num')

% % Display library function signatures
% libfunctionsview('QBladeDLL')

% % Call library function
% calllib('QBladeDLL',funcname,arguments)

%% Load OL data
clearvars;
% load('inputData\OLdata_steadyWind_irrWaves.mat')
load('inputData\OLdata_steadyWind_irrWaves_full.mat')

in1 = u;
out1 = y(:,1);
out2 = y(:,2);
dist1 = w;

%% Load wave data
% load('inputData\waveData.mat')
% waveData = mPitch(1:20:end);

%% DeePC parameters
N = 500; % lenght of data set
p = 20; % past data window
f = 20; % prediction window
Nbar = N-p-f+1;
i = 1;

controlParams.N = N;
controlParams.p = p;
controlParams.f = f;

%% Data detrending
mean_bladePitch = mean(in1);
uhat = in1 - mean_bladePitch;

mean_rotSpeed = mean(out1);
yhat1 = out1 - mean_rotSpeed;

mean_ptfmPitch = mean(out2);
yhat2 = out2 - mean_ptfmPitch;

mPitchHat = dist1; % already centered around 0

%% Scaling
% Scaling factors
uhat_max = 5; % Maximum expected input (deg)
% v0hat_max = 0.1*16; % Maximum expected wind disturbance (m/s)
mPitchHat_max = max(abs(mPitchHat)); % Maximum expected wave pitch moment (Nm)
omegaHat_max = 0.1* 12.1; % Maximum expected generator speed error (30% around steady-state value) (rpm)
pltfmPitchHat_max = 1.5; % Maximum expected platform pitch angle (deg)

% Signal scaling
u = uhat./uhat_max;
mPitch = mPitchHat./mPitchHat_max;
% waveData = waveData./mPitchHat_max;
y1 = yhat1./omegaHat_max;
y2 = yhat2./pltfmPitchHat_max;

y = [y1 y2];

%% Construct data matrices
% Use preview information?
% previewAns = input('Use preview information? y/n[y]: ','s');
% 
% while not(isempty(previewAns)) && not(strcmp(previewAns,'y')) && ...
%         not(strcmp(previewAns,'n'))    
%     disp('Invalid input.')
%     previewAns = input('Use preview information? y/n[y]: ','s');
% end
% 
% if isempty(previewAns) || strcmp(previewAns,'y')
%     previewFlag = 1;
% elseif strcmp(previewAns,'n')
%     previewFlag = 0;
% end
previewFlag = 1;
clear data;
% Past data
data.Up = constructHankelMat(u,i,p,Nbar);
data.Yp = constructHankelMat(y,i,p,Nbar);

% Future data
data.Uf = constructHankelMat(u,i+p,f,Nbar);
data.Yf = constructHankelMat(y,i+p,f,Nbar);

disturbMat = mPitch;

if previewFlag == 1
    data.Wp = constructHankelMat(disturbMat,i,p,Nbar); % past data
    data.Wf = constructHankelMat(disturbMat,i+p,f,Nbar); % future data
else
    data.Wp = [];
    data.Wf = [];
end

% Past data for prediction
data.uini = constructHankelMat(u,i+N-p,p,1);
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

%% Set up control loop
kFinal = 10000; % simulation steps
Ts = 0.05;
Ts_waves = 1;
Ts_ratio = Ts_waves/Ts;
kFinalWaves = kFinal/Ts_ratio;
tsim = 0:Ts:Ts*(kFinal-1);
nInputs = size(data.Up,1)/p;
nOutputs = size(data.Yp,1)/p;
if previewFlag == 1
    nDist = size(data.Wp,1)/p;
end

% Generate reference trajectory
ref_y1 = 12.12*ones(kFinal+f,1)';
ref_y1 = (ref_y1 - mean_rotSpeed)./omegaHat_max;
ref_y2 = zeros(kFinal+f,1)';
ref = [ref_y1; ref_y2];
ref = reshape(ref,[],1);

% Keep track of input and output sequences
uSeq = zeros(nInputs,kFinalWaves);
out = zeros(nOutputs,kFinalWaves);

%% Define control weights:

% weightOutputs diagonal matrix of size l-by-l, where l is the number of 
% output channels and the n-th element on the diagonal represents the 
% weight for the corresponding n-th output
weightOutputs = diag([1e2 0]);%preview
% weightOutputs =  1*diag(1);% no preview
controlParams.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of 
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
% weightInputs= 1e-3*diag(1); % no preview
weightInputs= 1e-1*diag(1); % preview
controlParams.R = kron(eye(f),weightInputs);

controlParams.lam_y = 1e3;
controlParams.lam_w = 1e3;

% Choose input bounds
bladePitch_min = 4;
bladePitch_max = 19;
controlParams.lbu = (bladePitch_min - mean_bladePitch)./uhat_max; % deg
controlParams.ubu = (bladePitch_max - mean_bladePitch)./uhat_max;

% Input rate constraint
duDeg = 8; % deg/s
duDeg = duDeg./uhat_max;
controlParams.duf = duDeg*Ts_waves;

% Output bounds
controlParams.lby = [-10; -0.75]; % rpm, degrees 
controlParams.lby = controlParams.lby./[omegaHat_max; pltfmPitchHat_max];
controlParams.uby = [10; 0.75];
controlParams.uby = controlParams.uby./[omegaHat_max; pltfmPitchHat_max];

% Measurement noise standard deviation
Std = 5e-5; 

% Choose optimization method
% method = input(['Optimization method: 1 - QP, ' '2 - SDP, ' '3 - NLP: ']);
method = 1; % QP

%% Control loop
calllib('QBladeDLL','createInstance',0,24)
calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')

% Offshore, steady wind, irregular waves
simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_steadyWind_IrrWaves.sim';

calllib('QBladeDLL','loadSimDefinition',simName)
calllib('QBladeDLL','initializeSimulation')

rotSpeedStr = 'Rotational Speed [rpm]';
mPitchStr = 'Y_g Mom. Diffraction [Nm]';
pitchStr = 'Pitch Angle Blade 1 [deg]';
ptfmPitchStr = 'NP Pitch Y_l [deg]';

ratedTq = 43093.55;

wt_out = zeros(nOutputs,kFinalWaves);
pitch_sim = zeros(kFinalWaves,3);
V_hub = zeros(kFinalWaves,3);
mPitchPreview = zeros(controlParams.f,1);

rotSpeedFull    = zeros(kFinal,1);
pitchFull       = zeros(kFinal,3);
timeFull        = zeros(kFinal,1);
ptfmPitchFull   = zeros(kFinal,1);
genPwrFull      = zeros(kFinal,1);
MytFull         = zeros(kFinal,1);
lssT            = zeros(kFinal,1);
Moop            = zeros(kFinal,1);

data.wf = zeros(controlParams.f,1);

% Turbulent wind:n
load('inputData\turbWind_16mps_long.mat')
v = turbWind;


%% Reach steady state
waitBar = waitbar(0,'Reaching steady-state') ;

for kSS = 1:1:1000
    calllib('QBladeDLL','advanceTurbineSimulation')

    % Set wind speed
    calllib('QBladeDLL','setPowerLawWind', v(kSS), 0, 0, 0, 87.6);

    theta_c = mean_bladePitch;
    calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);   

    waitbar(kSS/1000,waitBar,'Reaching steady-state')
end
close(waitBar)

%% Start buffer to get preview
waitBar = waitbar(0,'Preview buffer') ;

rot = zeros(controlParams.f*Ts_ratio,1);
bp = zeros(controlParams.f*Ts_ratio,1);

for kBuffer = 1:1:(controlParams.f*Ts_ratio)
    calllib('QBladeDLL','advanceTurbineSimulation')

    % Set wind speed
    calllib('QBladeDLL','setPowerLawWind', v(kSS+kBuffer), 0, 0, 0, 87.6);

    if mod(kBuffer-1, Ts_ratio) == 0
        idxPreview = floor(kBuffer/Ts_ratio) + 1;
        mPitchPreview_temp = calllib('QBladeDLL','getCustomData_at_num','Y_g Mom. Diffraction Offset [Nm]', 0, 0);
        mPitchPreview(idxPreview) = mPitchPreview_temp/mPitchHat_max;

        % Initialize trajectories for DeePC
        bp = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);
        uini_temp = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);
        uini_temp = (uini_temp - mean_bladePitch)/uhat_max;
        data.uini = [data.uini(nInputs+1:end); uini_temp]; 

        rot(kBuffer) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
        yini1_temp = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
        yini1_temp = (yini1_temp - mean_rotSpeed)./omegaHat_max;

        yini2_temp = calllib('QBladeDLL','getCustomData_at_num',ptfmPitchStr, 0, 0);
        yini2_temp = (yini2_temp - mean_ptfmPitch)./pltfmPitchHat_max;

        data.yini = [data.yini(nOutputs+1:end); yini1_temp; yini2_temp];

        if previewFlag == 1
            mPitch_temp = calllib('QBladeDLL','getCustomData_at_num',mPitchStr, 0, 0) ;
            mPitch_temp = mPitch_temp./mPitchHat_max;
            data.wini = [data.wini(nDist+1:end); mPitch_temp];
        end    
    end
    waitbar(kBuffer/((controlParams.f-1)*Ts_ratio),waitBar,'Getting preview data')
end
close(waitBar)

data.wf = mPitchPreview;

%% Start control
waitBar = waitbar(0,'Initializing Simulation') ;

for k = 1:1:kFinal
    fprintf("Iteration: %d \n", k)

    if mod(k-1,Ts_ratio) == 0  % once every second
        idxWaves = floor(k/Ts_ratio) + 1;        
        
        calllib('QBladeDLL','advanceTurbineSimulation')
         % Set wind speed
        calllib('QBladeDLL','setPowerLawWind', v(kSS+kBuffer+k), 0, 0, 0, 87.6);

        % Get rotor speed
        genPwrFull(k)     = calllib('QBladeDLL','getCustomData_at_num','Gen. Elec. Power [kW]', 0, 0);
        MytFull(k)        = calllib('QBladeDLL','getCustomData_at_num','Y_l Mom. TWR pos 0.000 [Nm]', 0, 0);
        lssT(k)           = calllib('QBladeDLL','getCustomData_at_num','Aero. LSS Torque [Nm]', 0, 0);
        Moop(k)           = calllib('QBladeDLL','getCustomData_at_num','Y_c RootBend. Mom. (OOP) BLD 1 [Nm]', 0, 0);

        wt_out(1,idxWaves) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
        rotSpeedFull(k) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
        wt_out(1,idxWaves) = wt_out(1,idxWaves) + Std.*randn(size(wt_out(1,idxWaves)));

        wt_out(2,idxWaves) = calllib('QBladeDLL','getCustomData_at_num',ptfmPitchStr, 0, 0);
        ptfmPitchFull(k) = calllib('QBladeDLL','getCustomData_at_num',ptfmPitchStr, 0, 0);
        wt_out(2,idxWaves) = wt_out(2,idxWaves) + Std.*randn(size(wt_out(2,idxWaves)));

        rotSpeed_out_deepc = (wt_out(1,idxWaves)-mean_rotSpeed)./omegaHat_max; % detrend and scale     
        ptfmPitch_out_deepc = (wt_out(2,idxWaves)-mean_ptfmPitch)./pltfmPitchHat_max; 
        data.yini = [data.yini(nOutputs+1:end); rotSpeed_out_deepc; ptfmPitch_out_deepc];

        % Reference trajectory
        rf = ref((k-1)*nOutputs+1:(k-1)*nOutputs+nOutputs*f);

        % Get preview
        if previewFlag == 1
            mPitchOffset = calllib('QBladeDLL','getCustomData_at_num','Y_g Mom. Diffraction Offset [Nm]', 0, 0);
            mPitchOffset = mPitchOffset/mPitchHat_max;
            data.wf = [data.wf(nDist+1:end); mPitchOffset];
        else
            data.wf = [];
        end

        % DeePC optimal control input        
        uStar = deepc2(data,rf,controlParams,method,ivFlag,previewFlag); 
        uSeq(:,idxWaves) = uStar*uhat_max + mean_bladePitch;

        % Set collective blade pitch
        theta_c = uSeq(:,idxWaves);
        calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);

        pitch_sim(idxWaves,:) = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);
        pitchFull(k,:)    = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);

        % Get wave forces
        mPitchHat_sim = calllib('QBladeDLL','getCustomData_at_num',mPitchStr, 0, 0) ;
        mPitch_sim = mPitchHat_sim./mPitchHat_max;

        data.uini = [data.uini(nInputs+1:end); uStar];        
        if previewFlag == 1
            data.wini = [data.wini(nDist+1:end); mPitch_sim];
        end        
    else
        calllib('QBladeDLL','advanceTurbineSimulation')
        % Set wind speed
        calllib('QBladeDLL','setPowerLawWind', v(kSS+kBuffer+k), 0, 0, 0, 87.6);

        calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);

        rotSpeedFull(k)   = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
        ptfmPitchFull(k)  = calllib('QBladeDLL','getCustomData_at_num',ptfmPitchStr, 0, 0);
        pitchFull(k,:)    = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);
        timeFull(k)       = calllib('QBladeDLL','getCustomData_at_num','Time [s]',0,0);
        genPwrFull(k)     = calllib('QBladeDLL','getCustomData_at_num','Gen. Elec. Power [kW]', 0, 0);
        MytFull(k)        = calllib('QBladeDLL','getCustomData_at_num','Y_l Mom. TWR pos 0.000 [Nm]', 0, 0);
        lssT(k)           = calllib('QBladeDLL','getCustomData_at_num','Aero. LSS Torque [Nm]', 0, 0);
        Moop(k)           = calllib('QBladeDLL','getCustomData_at_num','Y_c RootBend. Mom. (OOP) BLD 1 [Nm]', 0, 0);

    end

    waitbar(k/kFinal,waitBar,'Simulation Running')
end

close(waitBar)
calllib('QBladeDLL','closeInstance')

%% Plotting
% Controlled output
tsimWaves = 0:Ts_waves:(kFinalWaves-1)*Ts_waves;

figure
plot(tsimWaves,wt_out(1,:))
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
title('DeePC for reference tracking')
grid on
hold on
plot(tsimWaves,ref_y1(1:kFinalWaves)'*omegaHat_max+mean_rotSpeed,'k-') % reference
xline(Ts_waves*f,'k--','Future window size')
yline(controlParams.lby(1)*omegaHat_max + mean_rotSpeed,'r--')
yline(controlParams.uby(1)*omegaHat_max + mean_rotSpeed,'r--')
legend('Controlled output','Reference','Location','SouthEast')
set(gcf,'Color','White')

figure
plot(tsimWaves,wt_out(2,:))
xlabel('Time (in s)')
ylabel('Platform pitch (in deg)')
yline(controlParams.lby(2)*pltfmPitchHat_max + mean_ptfmPitch,'r--')
yline(controlParams.uby(2)*pltfmPitchHat_max + mean_ptfmPitch,'r--')
title('DeePC for reference tracking')
grid on
hold on
xline(Ts_waves*f,'k--','Future window size')
set(gcf,'Color','White')


% Control input sequence
figure
plot(tsimWaves,uSeq(1:length(tsimWaves)))
% hold on
% plot(tsimWaves,pitch_sim(:,3))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-5 25])
yline(controlParams.lbu*uhat_max + mean_bladePitch,'r--','LineWidth',1)
yline(controlParams.ubu*uhat_max + mean_bladePitch,'r--','LineWidth',1)
xline(Ts_waves*f,'k--','Future window size')
title('Control input')
grid on
% legend('commanded','actual')
set(gcf,'Color','White')

% Control input rate
figure
plot(tsimWaves(1:end-1),diff(uSeq(1:length(tsimWaves))*(1/Ts_waves)))
xlabel('Time (in s)')
ylabel('\theta_c rate (in deg/s)')
yline(8,'r--','LineWidth',1)
yline(-8,'r--','LineWidth',1)
xline(Ts_waves*f,'k--','Future window size')
ylim([-15 15])
title('Control input rate')
grid on
set(gcf,'Color','White')

%% Save
% save('outputData\turb_np_noss.mat','tsim','rotSpeed','ref','kFinal','Ts','f','uSeq','controlParams')
% save('outputData\turb_np_noss_2.mat','tsim','rotSpeed','ref','kFinal','Ts','f','uSeq','controlParams')
% save('outputData\turb_wpStatic_noss_1.mat','tsim','rotSpeed','ref','kFinal','Ts','f','uSeq','controlParams')
% save('outputData\turb_wpDynamic_noss_1.mat','tsim','rotSpeed','ref','kFinal','Ts','f','uSeq','controlParams')
% 
% h = findobj('type','figure');
% savefig(flip(h),'figures\wind_preview.fig')

% save('outputData\previewDeePC.mat','rotSpeedFull','timeFull','pitchFull');
save('outputDataFinal\turb_deepc1.mat')
