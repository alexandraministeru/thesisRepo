%  Simulation of NREL5MW+OC3SPAR in steady wind 16m/s and still water.
%  The purpose is to reach a steady state and extract the value of the
%  blade pitch command to later use during OL data  collection for DeePC.

%% Clean environment
clearvars;close all; clc
rng('default')
cd 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\MATLAB_Files\waveFF\'
addpath(genpath('functions'))

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
%
% % Display library function signatures
% libfunctionsview('QBladeDLL')
%
% % Call library function
% calllib('QBladeDLL',funcname,arguments)

%% OL simulation
Ts_waves = 1;
Ts_wt = 0.05;
Ts_ratio = Ts_waves/Ts_wt;

calllib('QBladeDLL','createInstance',0,24)
calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')
 
% Offshore, steady wind, still water
simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_SPAR_irrWaves_steadyWind.sim';
calllib('QBladeDLL','loadSimDefinition',simName)

% % Set timestep size
% Ts = 0.5;
% calllib('QBladeDLL','setTimestepSize', Ts);

% Initialize simulation
calllib('QBladeDLL','initializeSimulation')

simTime = 700*Ts_ratio; %in timestep, actual time is timestep*#timesteps

% Generate input data: wind speed and collective blade pitch
simTimeWaves = simTime/Ts_ratio;
tWaves = 0:Ts_waves:(simTimeWaves-1)*Ts_waves; % time vector
tWind = 0:Ts_wt:(simTime-1)*Ts_wt;

u_bladePitch = idinput(length(tWaves),'PRBS',[0 1/10],[-2 2]); % in degrees
u_bladePitch = u_bladePitch + 11.6679; 

% % Generate horizontal wind speed disturbance input
u_windSpeed = 16 + 0.5.*randn(length(tWind),1);

ratedTq = 43093.55; % Nm

rotSpeedStr = 'Rotational Speed [rpm]';
pitchStr = 'Pitch Angle Blade 1 [deg]';
mPitchStr = 'Y_g Mom. Diffraction [Nm]';
ptfmPitchStr = 'NP Pitch Y_l [deg]';

%% Reach steady-state
kInit = 2000; % 100 seconds
mPitch_check = zeros(kInit,1);
pitchCheck = zeros(kInit,1);
V_hub_check = zeros(kInit,3);
rotSpeed_check = zeros(kInit,1);
ptfmPitch_check = zeros(kInit,1);
waitBar = waitbar(0,'Reaching steady-state') ;

% Reach a steady state first
for idxInit = 1:1:kInit
    calllib('QBladeDLL','advanceTurbineSimulation')

    % % Read pitch angle and horizontal wind speed at hub
    mPitch_check(idxInit) = calllib('QBladeDLL','getCustomData_at_num',mPitchStr, 0, 0) ;
    ptfmPitch_check(idxInit) = calllib('QBladeDLL','getCustomData_at_num',ptfmPitchStr, 0, 0) ;
    pitchCheck(idxInit,:) = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);  
    V_hub_check(idxInit,:) = calllib('QBladeDLL','getWindspeed', 0, 0, 87.6, [0 0 0]);
    rotSpeed_check(idxInit,:) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0) ;    

    % Set steady-state value of collective blade pitch angle
    theta_c = 11.6679;
    calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);      
    
    waitbar(idxInit/kInit,waitBar,'Reaching steady-state')
end
close(waitBar)

tInit = 0:Ts_wt:(kInit-1)*Ts_wt;

figure()
plot(tInit,pitchCheck,'LineWidth',1)
grid on
xlabel('Time (in s)')
ylabel('Pitch angle (in deg)')
title('Collective blade pitch')
grid on
set(gcf,'Color','White')

figure()
plot(tInit,ptfmPitch_check,'LineWidth',1)
grid on
xlabel('Time (in s)')
ylabel('Platform pitch angle (in deg)')
title('Platform pitch angle')
grid on
set(gcf,'Color','White')

figure()
plot(tInit,mPitch_check,'LineWidth',1)
grid on
xlabel('Time (in s)')
ylabel('Pitch moment (in Nm)')
title('Platform pitch moment')
grid on
set(gcf,'Color','White')

figure()
plot(tInit,V_hub_check(:,1),'LineWidth',1)
grid on
xlabel('Time (in s)')
ylabel('Horizontal wind speed (in m/s)')
title('Wind speed at hub')
% legend('WindVelX','WindVelY','WindVelY',Location='southeast')
grid on
set(gcf,'Color','White')

figure()
plot(tInit,rotSpeed_check,'LineWidth',1)
grid on
xlabel('Time (in s)')
ylabel('Wind speed (in rpm)')
title('Rotor speed')
grid on
set(gcf,'Color','White')

%% Collect OL data
kFinalWaves = simTime/Ts_ratio;
rotSpeed = zeros(kFinalWaves,1);
mPitch = zeros(kFinalWaves,1);
ptfmPitch = zeros(kFinalWaves,1);
pitch = zeros(kFinalWaves,1);
pitchFull = zeros(simTime,1);
V_hub = zeros(kFinalWaves,3);
rotSpeedFull = zeros(simTime,1);
ptfmPitchFull = zeros(simTime,1);
waitBar = waitbar(0,'Initializing Simulation') ;

idxWaves = 0;

for idxWT = 1:1:simTime   
    calllib('QBladeDLL','advanceTurbineSimulation')

    % % Set wind speed
    % calllib('QBladeDLL','setPowerLawWind', u_windSpeed(idxWT), 0, 0, 0, 87.6);

    waveFlag = mod(idxWT-1,Ts_ratio);    
    
    if waveFlag == 0
        idxWaves = idxWaves + 1;        
    
        % Get wave forces
        mPitch(idxWaves) = calllib('QBladeDLL','getCustomData_at_num',mPitchStr, 0, 0) ;
    
        % Get rotor speed and pitch angle
        rotSpeed(idxWaves,:) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0) ;   
        ptfmPitch(idxWaves,:) = calllib('QBladeDLL','getCustomData_at_num',ptfmPitchStr, 0, 0) ;
        pitch(idxWaves,:) = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);  
        V_hub(idxWaves,:) = calllib('QBladeDLL','getWindspeed', 0, 0, 87.6, [0 0 0]);

        % Set collective blade pitch
        theta_c = u_bladePitch(idxWaves);
        calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);

    end
    rotSpeedFull(idxWT,:) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0) ;
    pitchFull(idxWT,:) = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);  
    ptfmPitchFull(idxWT,:) = calllib('QBladeDLL','getCustomData_at_num',ptfmPitchStr, 0, 0);  
    waitbar(idxWT/simTime,waitBar,'Simulation Running')
end
close(waitBar)
calllib('QBladeDLL','closeInstance')

% Add measurement noise
% Std = 5e-5; % measurement noise standard deviation
Std = 1e-4; % measurement noise standard deviation
rotSpeed = rotSpeed + Std.*randn(size(rotSpeed));

%%
figure
plot(tWaves,rotSpeed)
grid on
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
xlim('tight')
legend('Rotor speed','Location','southeast')
set(gcf,'Color','White')
hold on

figure
plot(tWaves,u_bladePitch)
hold on
plot(tWaves,pitch)
grid on
xlabel('Time (in s)')
ylabel(pitchStr)
xlim('tight')
legend('Commanded','Actual subsampled')
set(gcf,'Color','White')

figure
plot(tWaves,mPitch)
grid on
xlabel('Time (in s)')
ylabel(mPitchStr)
set(gcf,'Color','White')

figure
plot(tWaves,ptfmPitch)
grid on
xlabel('Time (in s)')
ylabel('Horizontal wind speed')
set(gcf,'Color','White')

figure
plot(tWaves,V_hub(:,1))
grid on
xlim('tight')
xlabel('Time (in s)')
ylabel('Horizontal wind speed (in m/s)')
legend('Set wind speed at hub','Actual wind speed at hub','Location','southeast')
set(gcf,'Color','White')

%%
figure
plot(0:Ts_wt:(simTime-1)*Ts_wt,rotSpeedFull)
grid on
xlabel('Time (in s)')
ylabel('Rotor speed real samples(in rpm)')
xlim('tight')
legend('Rotor speed','Location','southeast')
set(gcf,'Color','White')
hold on

figure
plot(0:Ts_wt:(simTime-1)*Ts_wt,pitchFull)
grid on
xlabel('Time (in s)')
ylabel('Blade pich angle real samples(in deg)')
xlim('tight')
legend('Pitch angle','Location','southeast')
set(gcf,'Color','White')
hold on

figure
plot(0:Ts_wt:(simTime-1)*Ts_wt,ptfmPitchFull)
grid on
xlabel('Time (in s)')
ylabel('Platform pich angle real samples(in deg)')
xlim('tight')
legend('Pitch angle','Location','southeast')
set(gcf,'Color','White')
hold on

%%
u = pitch; % (100:end);
y = [rotSpeed ptfmPitch]; % (100:end);
% w = [V_hub(:,1) mPitch]; % (100:end);

w = mPitch; % (100:end);

save('inputData\OLdata_steadyWind_irrWaves_full_1e-4.mat','u','y','w')

% 
% save('inputData\OLdata_steadyWind_irrWaves_fullSys.mat','tWind','tWaves',...
%     'mPitch_check','rotSpeed_check','rotSpeedFull','pitchFull',...
%     'u_bladePitch','Std')