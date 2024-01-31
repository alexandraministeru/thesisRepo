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
% QBladeFunctionHelp('advanceController_at_num')

% % Display library function signatures
% libfunctionsview('QBladeDLL')

% % Call library function
% calllib('QBladeDLL',funcname,arguments)

%% OL simulation
calllib('QBladeDLL','createInstance',0,24)
calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')

% From a project file:
% % projectFile = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL5MW.qpr';
% % calllib('QBladeDLL','loadProject',projectFile)
 
% % Onshore, turbulent wind
% simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL5MW\NREL_5MW_turbWind.sim';

% Offshore, turbulent wind, still water
simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_SPAR_stillWater_turbWind.sim';

% % Offshore, turbulent wind, irregular waves
% simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_SPAR_irrWaves_turbWind.sim';

% Offshore, steady wind, irregular waves
% simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_SPAR_irrWaves_steadyWind.sim';

calllib('QBladeDLL','loadSimDefinition',simName)
calllib('QBladeDLL','initializeSimulation')

Ts = 0.05;
simTime = 2000; %in timestep, actual time is timestep*#timesteps
rotSpeedStr = 'Rotational Speed [rpm]';
pitchStr = 'Pitch Angle Blade 1 [deg]';

% Generate input data: wind speed and collective blade pitch
t = 0:Ts:(simTime-1)*Ts;

u_bladePitch = idinput(length(t),'PRBS',[0 1/10],[-2 2]); % in degrees
u_bladePitch = u_bladePitch + 11.79; 

% % Generate horizontal wind speed disturbance input
u_windSpeed = 16 + 0.5.*randn(length(t),1);

ratedTq = 43093.55;

% Reach a steady state first
% for kInit = 1:1500
%     calllib('QBladeDLL','advanceTurbineSimulation')
%     theta_c = 11.79;
%     calllib('QBladeDLL','setPowerLawWind', 16, 0, 0, 0, 87.6); % Steady wind
%     calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);    
% end

rotSpeed = zeros(simTime,1);
pitch = zeros(simTime,1);
V_hub = zeros(simTime,3);
waitBar = waitbar(0,'Initializing Simulation') ;

% calllib('QBladeDLL','advanceTurbineSimulation')

for i = 1:1:simTime
    % Set wind speed
    calllib('QBladeDLL','setPowerLawWind', u_windSpeed(i), 0, 0, 0, 90);

    % Set collective blade pitch
    theta_c = u_bladePitch(i);
    calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);

    calllib('QBladeDLL','advanceTurbineSimulation')

    % Get rotor speed
    rotSpeed(i,:) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0) ;
    pitch(i,:) = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);    

    V_hub(i,:) = calllib('QBladeDLL','getWindspeed', 0, 0, 90, [0 0 0]);
    waitbar(i/simTime,waitBar,'Simulation Running')
end
close(waitBar)
calllib('QBladeDLL','closeInstance')

% Add measurement noise
std = 5e-5; % measurement noise standard deviation
rotSpeed = rotSpeed + std.*randn(size(rotSpeed));

figure
plot(t,rotSpeed)
grid on
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
legend('Rotor speed','Location','southeast')
set(gcf,'Color','White')
hold on

figure
plot(t,pitch)
grid on
xlabel('Time (in s)')
ylabel(pitchStr)
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

%% Save OL data
% save('inputData\openLoopData_sim_noSS.mat','rotSpeed','u_bladePitch','u_windSpeed','t')
% save('inputData\openLoopData_sim_SS.mat','rotSpeed','u_bladePitch','u_windSpeed','t')
save('inputData\openLoopData_oc3spar_stillWater.mat','rotSpeed','u_bladePitch','u_windSpeed','t')
% save('inputData\openLoopData_oc3spar_irrWaves.mat','rotSpeed','u_bladePitch','u_windSpeed','t')
% load('inputData\openLoopData.mat','rotSpeed','u_bladePitch','u_windSpeed')
% load('inputData\openLoopData_oc3spar.mat','rotSpeed','u_bladePitch','u_windSpeed')

%% Load OL data
clearvars;
% load('inputData\openLoopData_sim.mat','rotSpeed','u_bladePitch','u_windSpeed')
% load('inputData\openLoopData_sim_noSS.mat','rotSpeed','u_bladePitch','u_windSpeed')
load('inputData\openLoopData_oc3spar_stillWater.mat','rotSpeed','u_bladePitch','u_windSpeed')
% load('inputData\openLoopData_oc3spar_irrWaves.mat','rotSpeed','u_bladePitch','u_windSpeed')


%% DeePC parameters
N = 500; % lenght of data set
p = 40; % past data window
f = 20; % prediction window
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

y = rotSpeed;
clear data;
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

% Past data for prediction
data.uini = constructHankelMat(u_bladePitch,i+N-p,p,1);
data.yini = constructHankelMat(y,i+N-p,p,1);
if previewFlag == 1
    data.wini = constructHankelMat(u_windSpeed,i+N-p,p,1);
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
kFinal = 600; % simulation steps
Ts = 0.05;
tsim = 0:Ts:Ts*(kFinal-1);
nInputs = size(data.Up,1)/p;
nOutputs = size(data.Yp,1)/p;
if previewFlag == 1
    nDist = size(data.Wp,1)/p;
end

% Generate reference trajectory
ref = 12.1*ones(kFinal+f,1);
% ref(200:end) = 100; % step in reference

% Keep track of input and output sequences
uSeq = zeros(nInputs,kFinal);
out = zeros(nOutputs,kFinal);

% Define control weights:
%
% weightOutputs diagonal matrix of size l-by-l, where l is the number of 
% output channels and the n-th element on the diagonal represents the 
% weight for the corresponding n-th output

if ivFlag == 1    
    % weightOutputs = 1e5*diag(1);%preview
    % weightOutputs = 1e3*diag(1);% no preview
    weightOutputs = 1*diag(1);% no preview
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
controlParams.lbu = 0; % deg
controlParams.ubu = 22;

% Input rate constraint
duDeg = 8; % deg/s
controlParams.duf = duDeg*Ts;

% Measurement noise standard deviation
std = 5e-5;

% Choose optimization method
% method = input(['Optimization method: 1 - QP, ' '2 - SDP, ' '3 - NLP: ']);
method = 1;

% Wind disturbance
% Extreme operating gust
% load('inputData\eog_16mps.mat','Wind1VelX','Time')
% v = interp1(Time,Wind1VelX,tsim)'; % resample with current sampling period
% v = [v;16*ones(f,1)];

% Turbulent wind
% load('inputData\turbWind_16mps_long.mat') %turbulent wind obtained from a previous FAST simulation
% v = turbWind;

% % Turbulent wind2
% load('outputData\turbWind2.mat') %turbulent wind obtained from a previous QBlade simulation
% v = V_hub(:,1);

% load('inputData\turbWind_QBlade_onshore.mat')
% v = windData;

%% Test wind preview
% calllib('QBladeDLL','createInstance',0,24)
% calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')
% 
% simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL5MW\NREL_5MW_turbWind.sim';
% calllib('QBladeDLL','loadSimDefinition',simName)
% 
% calllib('QBladeDLL','initializeSimulation')
% 
% V_hub_pre = zeros(10,3);
% V_hub_post = zeros(10,3); 
% calllib('QBladeDLL','advanceTurbineSimulation')
% 
% for k=1:10
%     V_hub_pre(k,:) = calllib('QBladeDLL','getWindspeed', 0, 0, 87.6, [0 0 0]);    
%     calllib('QBladeDLL','advanceTurbineSimulation')
%     V_hub_post(k,:) = calllib('QBladeDLL','getWindspeed', 0, 0, 87.6, [0 0 0]);
% end
% 
% calllib('QBladeDLL','closeInstance')
% 
% figure
% plot(V_hub_post(:,1),'*-')
% hold on
% plot(V_hub_post2(:,1),'o-')
% legend('pre','post')

%% Control loop
calllib('QBladeDLL','createInstance',0,24)
calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')

% From a project file:
% % projectFile = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL5MW.qpr';
% % calllib('QBladeDLL','loadProject',projectFile)
 
% % % Onshore, turbulent wind
% simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL5MW\NREL_5MW_turbWind.sim';

% Offshore, turbulent wind, still water
simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_SPAR_stillWater_turbWind.sim';

% % Offshore, turbulent wind, irregular waves
% simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_SPAR_irrWaves_turbWind.sim';

% % Offshore, steady wind, irregular waves
% simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_SPAR_irrWaves_steadyWind.sim';

calllib('QBladeDLL','loadSimDefinition',simName)
calllib('QBladeDLL','initializeSimulation')

rotSpeedStr = 'Rotational Speed [rpm]';
pitchStr = 'Pitch Angle Blade 1 [deg]';

% pitchStr = 'Pitch Angle Blade 1 [deg]';
% for kInit = 1:1000
%     calllib('QBladeDLL','setPowerLawWind', 16, 0, 0, 0, 87.6);
%     % Set collective blade pitch
%     theta_c = 11.79;
%     calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);
% 
%     calllib('QBladeDLL','advanceTurbineSimulation')
% 
%     rotSpeedInit(:,kInit) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
%     pitchInit(:,kInit) = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 1, 0);
% end
% calllib('QBladeDLL','closeInstance')
% figure
% plot(rotSpeedInit)
% figure
% plot(pitchInit)

waitBar = waitbar(0,'Initializing Simulation') ;
rotSpeed_out = zeros(nOutputs,kFinal);
pitch = zeros(1,kFinal);
ratedTq = 43093.55;
control_out = zeros(5,kFinal);
V_hub = zeros(kFinal,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Preview:
V_hub_preview = zeros(f,3);
avgWindSpeed = 16; % m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calllib('QBladeDLL','advanceTurbineSimulation')

for k = 1:1:kFinal
    % % "Start-up"
    % if k == 1
    %     for kInit = 1:500
    %         calllib('QBladeDLL','setPowerLawWind', 16, 0, 0, 0, 87.6);
    %         % Set collective blade pitch
    %         theta_c = 11.79;
    %         calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);
    %         calllib('QBladeDLL','advanceTurbineSimulation')          
    %     end
    % end

    fprintf("Iteration: %d \n", k)

    % % Set wind speed
    % calllib('QBladeDLL','setPowerLawWind', v(k), 0, 0, 0, 87.6);    

    % Reference trajectory
    rf = ref((k-1)*nOutputs+1:(k-1)*nOutputs+nOutputs*f);

    % Wind preview
    if previewFlag == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Get preview
        for idxPreview = 0:1:f-1
            xDist = idxPreview*Ts*avgWindSpeed;
            V_hub_preview(idxPreview+1,:) = calllib('QBladeDLL','getWindspeed', -xDist, 0, 90, [0 0 0]);
        end
        data.wf = V_hub_preview(:,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % data.wf = v(k:k+f-1); 
    else
        data.wf = [];
    end

    % DeePC optimal control input
    uStar = deepc(data,rf,controlParams,method,ivFlag,previewFlag);
    uSeq(:,k) = uStar;   

    % Set collective blade pitch
    theta_c = uSeq(:,k);
    calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);

    % Get horizontal wind speed at hub
    V_hub(k,:) = calllib('QBladeDLL','getWindspeed', 0, 0, 90, [0 0 0]);

    calllib('QBladeDLL','advanceTurbineSimulation')
    control_out(:,k) = calllib('QBladeDLL','advanceController_at_num',[0 0 0 0 0], 0);

    % Get rotor speed
    pitch(:,k)    = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 1, 0);
    rotSpeed_out(:,k) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
    rotSpeed_out(:,k) = rotSpeed_out(:,k) + std.*randn(size(rotSpeed_out(:,k)));
    
          
   % Update past data with most recent I/O data
    data.uini = [data.uini(nInputs+1:end); uSeq(:,k)];
    data.yini = [data.yini(nOutputs+1:end); rotSpeed_out(:,k)];
    if previewFlag == 1
        % data.wini = [data.wini(nDist+1:end); v(k)];
        data.wini = [data.wini(nDist+1:end); V_hub(k,1)];
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
plot(tsim,ref(1:kFinal)) % reference
xline(Ts*f,'k--','Future window size')
% xline(Ts*stepIdxs(1),'k--','Reference step')
legend('Controlled output','Reference','Location','SouthEast')
set(gcf,'Color','White')

figure
plot(tsim,pitch)
xlabel('Time (in s)')
ylabel('Pitch angle')
title('Pitch angle')
grid on
hold on
xline(Ts*f,'k--','Future window size')
legend('Controlled output','Reference','Location','SouthEast')
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

% Wind disturbance
figure
plot(tsim,v(1:length(tsim)))
hold on
plot(tsim,V_hub(:,1))
xlabel('Time (in s)')
ylabel('Wind speed (in m/s)')
title('Horizontal wind speed around 16mps linearization point')
grid on
set(gcf,'Color','White')

%% Save
% save('outputData\turb_np_noss.mat','tsim','rotSpeed','ref','kFinal','Ts','f','uSeq','controlParams')
% save('outputData\turb_np_noss_2.mat','tsim','rotSpeed','ref','kFinal','Ts','f','uSeq','controlParams')
% save('outputData\turb_wpStatic_noss_1.mat','tsim','rotSpeed','ref','kFinal','Ts','f','uSeq','controlParams')
% save('outputData\turb_wpDynamic_noss_1.mat','tsim','rotSpeed','ref','kFinal','Ts','f','uSeq','controlParams')


