%%
% DeePC FF control

%% Clean environment
clearvars;close all; clc
rng('default')
cd 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\MATLAB_Files\ffControl\'

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
% QBladeFunctionHelp('getCustomData_at_num')

% % Display library function signatures
% libfunctionsview('QBladeDLL')

% % Call library function
% calllib('QBladeDLL',funcname,arguments)

%% OL simulation
calllib('QBladeDLL','createInstance',0,24)
calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')
 
% Offshore, steady wind, irregular waves
simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_steadyWind_IrrWaves.sim';

calllib('QBladeDLL','loadSimDefinition',simName)
calllib('QBladeDLL','initializeSimulation')

Ts = 0.05;
simTime = 16000; %in timestep, actual time is timestep*#timesteps
rotSpeedStr = 'Rotational Speed [rpm]';
mPitchStr = 'Y_g Mom. Diffraction [Nm]';
pitchStr = 'Pitch Angle Blade 1 [deg]';

% Generate input data: wind speed and collective blade pitch
t = 0:Ts:(simTime-1)*Ts;

u_bladePitch = idinput(length(t),'PRBS',[0 1/200],[-2 2]); % in degrees
u_bladePitch = u_bladePitch + 11.79; 

ratedTq = 43093.55;

% Reach a steady state first
% for kInit = 1:1500
%     calllib('QBladeDLL','advanceTurbineSimulation')
%     theta_c = 11.79;
%     calllib('QBladeDLL','setPowerLawWind', 16, 0, 0, 0, 87.6); % Steady wind
%     calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);    
% end

rotSpeed = zeros(simTime,1);
mPitch = zeros(simTime,1);
pitch = zeros(simTime,1);
V_hub = zeros(simTime,3);
waitBar = waitbar(0,'Initializing Simulation') ;

for i = 1:1:simTime   
    % Set collective blade pitch
    theta_c = u_bladePitch(i);
    calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);
    % Get wave forces
    mPitch(i) = calllib('QBladeDLL','getCustomData_at_num',mPitchStr, 0, 0) ;

    calllib('QBladeDLL','advanceTurbineSimulation')

    % Get rotor speed
    rotSpeed(i,:) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0) ;    
    pitch(i,:) = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);  

    waitbar(i/simTime,waitBar,'Simulation Running')
end
close(waitBar)
calllib('QBladeDLL','closeInstance')

% Add measurement noise
Std = 5e-3; % measurement noise standard deviation
rotSpeed = rotSpeed + Std.*randn(size(rotSpeed));

figure
plot(t,rotSpeed)
grid on
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
legend('Rotor speed','Location','southeast')
set(gcf,'Color','White')
hold on

figure
plot(t,u_bladePitch)
grid on
xlabel('Time (in s)')
ylabel(pitchStr)
set(gcf,'Color','White')

figure
plot(t,mPitch)
grid on
xlabel('Time (in s)')
ylabel(mPitchStr)
set(gcf,'Color','White')

figure
plot(t,u_windSpeed)
hold on
plot(t,V_hub(:,1))
grid on
xlabel('Time (in s)')
ylabel('Horizontal wind speed (in m/s)')
legend('Wind Speed at hub','Location','southeast')
set(gcf,'Color','White')

Ts_waves = 1;
Ts_ratio = Ts_waves/Ts;

rotSpeed = rotSpeed(1:Ts_ratio:end);
u_bladePitch = u_bladePitch(1:Ts_ratio:end);
mPitch = mPitch(1:Ts_ratio:end);
t = t(1:Ts_ratio:end);

%% Save OL data
save('inputData\openLoopData_waves_ts1.mat','rotSpeed','u_bladePitch','mPitch','t')

%% Load OL data
clearvars;
load('inputData\openLoopData_waves_ts1.mat','rotSpeed','u_bladePitch','mPitch')

%% DeePC parameters
N = 600; % lenght of data set
p = 20; % past data window
f = 20; % prediction window
Nbar = N-p-f+1;
i = 1;

controlParams.N = N;
controlParams.p = p;
controlParams.f = f;

%% Scaling
% % Scaling factors
% uhat_max = 10; % Maximum expected input (deg)
% v0hat_max = 0.1*16; % Maximum expected wind disturbance (m/s)
% mPitchHat_max = max(mPitch); % Maximum expected wave pitch moment (Nm)
% 
% % ehat_max = 0.1* 1173.69; % Maximum expected generator speed error (10% around linearization OP) (rpm)
% omegaHat_max = 0.1* 12.1; % Maximum expected generator speed error (10% around linearization OP) (rpm)
% 
% % Scale signals
% uhat = u_bladePitch - mean(u_bladePitch);
% u = uhat./uhat_max;
% 
% mPitchHat = mPitch;
% mPitch = mPitchHat./mPitchHat_max;
% 
% yhat = rotSpeed - mean(rotSpeed);
% y = yhat./omegaHat_max;

u = u_bladePitch;
y = rotSpeed;
%% Construct data matrices
% Use preview information?
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
kFinal = 2000; % simulation steps
Ts = 0.05;
tsim = 0:Ts:Ts*(kFinal-1);
nInputs = size(data.Up,1)/p;
nOutputs = size(data.Yp,1)/p;
if previewFlag == 1
    nDist = size(data.Wp,1)/p;
end

% Generate reference trajectory
ref = 12.1*ones(kFinal+f,1);
% ref = (ref-mean(rotSpeed))./omegaHat_max;

% Keep track of input and output sequences
uSeq = zeros(nInputs,kFinal);
out = zeros(nOutputs,kFinal);

% Define control weights:
%
% weightOutputs diagonal matrix of size l-by-l, where l is the number of 
% output channels and the n-th element on the diagonal represents the 
% weight for the corresponding n-th output


% weightOutputs = 1e5*diag(1);%preview
weightOutputs = 1e3*diag(1);% no preview
controlParams.Q = kron(eye(f),weightOutputs);

% weightInputs diagonal matrix of size m-by-m, where m is the number of 
% input channels and the n-th element on the diagonal represents the weight
% for the corresponding n-th input
weightInputs= diag(1); 
controlParams.R = kron(eye(f),weightInputs);

% Choose input bounds
% controlParams.lbu = (0-mean(u_bladePitch))./uhat_max; % deg
% controlParams.ubu = (22-mean(u_bladePitch))./uhat_max;
controlParams.lbu = 0; % deg
controlParams.ubu = 22;

% Input rate constraint
duDeg = 8; % deg/s
% duDeg = duDeg./uhat_max;
controlParams.duf = duDeg*1;

% Measurement noise standard deviation
Std = 5e-3; 

% Choose optimization method
% method = input(['Optimization method: 1 - QP, ' '2 - SDP, ' '3 - NLP: ']);
method = 1;

%% Control loop
calllib('QBladeDLL','createInstance',0,24)
calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')

% Offshore, steady wind, irregular waves
simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_steadyWind_IrrWaves.sim';

calllib('QBladeDLL','loadSimDefinition',simName)
calllib('QBladeDLL','initializeSimulation')

rotSpeedStr = 'Rotational Speed [rpm]';
mPitchStr = 'Y_g Mom. Diffraction [Nm]';

waitBar = waitbar(0,'Initializing Simulation') ;
rotSpeed_out = zeros(nOutputs,kFinal);
ratedTq = 43093.55;
control_out = zeros(5,kFinal);
V_hub = zeros(kFinal,3);

mPitchPreview = zeros(controlParams.f,1);
Ts_ratio = 20;

calllib('QBladeDLL','advanceTurbineSimulation')

for k = 1:1:kFinal
    if mod(k-1,Ts_ratio) == 0 && k > controlParams.f

        fprintf("Iteration: %d \n", k)

        % Reference trajectory
        rf = ref((k-1)*nOutputs+1:(k-1)*nOutputs+nOutputs*f);

        % Get preview
        if previewFlag == 1
            mPitchOffset = calllib('QBladeDLL','getCustomData_at_num','Y_g Mom. Diffraction Offset [Nm]', 0, 0);
            % mPitchOffset = mPitchOffset/mPitchHat_max;
            data.wf = [mPitchPreview(2:end); mPitchOffset];
        else
            data.wf = [];
        end

        % DeePC optimal control input
        uStar = deepc(data,rf,controlParams,method,ivFlag,previewFlag);
        uSeq(:,k) = uStar;%*uhat_max + mean(u_bladePitch);

        % Set collective blade pitch
        theta_c = uSeq(:,k);
        calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);

        % Get wave forces
        mPitchHat_sim = calllib('QBladeDLL','getCustomData_at_num',mPitchStr, 0, 0) ;
        mPitch_sim = mPitchHat_sim;%./mPitchHat_max;

        calllib('QBladeDLL','advanceTurbineSimulation')
        control_out(:,k) = calllib('QBladeDLL','advanceController_at_num',[0 0 0 0 0], 0);

        % Get rotor speed
        rotSpeed_out(:,k) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
        rotSpeed_out(:,k) = rotSpeed_out(:,k) + Std.*randn(size(rotSpeed_out(:,k)));

        % rotSpeed_out_deepc = (rotSpeed_out(:,k)-mean(rotSpeed))./omegaHat_max;

        % Update past data with most recent I/O data
        data.uini = [data.uini(nInputs+1:end); uStar];
        data.yini = [data.yini(nOutputs+1:end); rotSpeed_out(:,k)];
        if previewFlag == 1
            data.wini = [data.wini(nDist+1:end); mPitch_sim];
        end

        waitbar(k/kFinal,waitBar,'Simulation Running')
    elseif k <= controlParams.f
        mPitchPreview_temp = calllib('QBladeDLL','getCustomData_at_num','Y_g Mom. Diffraction Offset [Nm]', 0, 0);
        mPitchPreview(k) = mPitchPreview_temp;%/mPitchHat_max;
    end
    waitbar(k/kFinal,waitBar,'Simulation Running')
end

close(waitBar)
calllib('QBladeDLL','closeInstance')

%% Plotting
% Controlled output
figure
plot(tsim,rotSpeed_out)
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
title('DeePC for reference tracking')
grid on
hold on
% plot(tsim,ref(1:kFinal)*omegaHat_max+mean(rotSpeed)) % reference
xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
% legend('Controlled output','Reference','Location','SouthEast')
set(gcf,'Color','White')

% Control input sequence
figure
plot(tsim,uSeq)
hold on
plot(tsim,control_out(3,:))
xlabel('Time (in s)')
ylabel('\theta_c (in deg)')
ylim([-5 25])
yline(controlParams.lbu,'r--','LineWidth',1)
yline(controlParams.ubu,'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input')
grid on
set(gcf,'Color','White')

% Control input rate
figure
plot(tsim(1:end-1),diff(uSeq)*(1/Ts))
xlabel('Time (in s)')
ylabel('\theta_c rate (in deg/s)')
yline(duDeg,'r--','LineWidth',1)
yline(-duDeg,'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
ylim([-15 15])
title('Control input rate')
grid on
set(gcf,'Color','White')


