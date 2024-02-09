%% Clean environment
clearvars;close all; clc
rng('default')
cd 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\MATLAB_Files\waveFF\'

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

%%
calllib('QBladeDLL','createInstance',0,24)
calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')

simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_steadyWind_IrrWaves_pi1.sim';
calllib('QBladeDLL','loadSimDefinition',simName)
calllib('QBladeDLL','initializeSimulation')
Ts = 0.05;
simTime = 15600; % steady
% simTime = 10000; % turbulent
t = 0:Ts:(simTime-1)*Ts;

pitchFull       = zeros(simTime,3);
rotSpeedFull    = zeros(simTime,1);
ptfmPitchFull   = zeros(simTime,1);
genPwrFull      = zeros(simTime,1);
MytFull         = zeros(simTime,1);
lssT            = zeros(simTime,1);
timeFull        = zeros(simTime,1);
Moop            = zeros(simTime,1);

rotSpeedStr = 'Rotational Speed [rpm]';
pitchStr = 'Pitch Angle Blade 1 [deg]';
ptfmPitchStr = 'NP Pitch Y_l [deg]';
% 'Gen. Elec. Power [kW]'
% 'Y_l Mom. TWR pos 0.000 [Nm]'
% 'Aero. LSS Torque [Nm]'
% 'Y_c RootBend. Mom. (OOP) BLD 1 [Nm]'

load('inputData\turbWind_16mps_long.mat');
v = turbWind;

waitBar = waitbar(0,'Preview buffer') ;

for idx = 1:1:simTime
    disp(idx)
    calllib('QBladeDLL','advanceTurbineSimulation')
    calllib('QBladeDLL','advanceController_at_num', ones(1,5), 0)   

    %  % Set wind speed
    % calllib('QBladeDLL','setPowerLawWind', v(idx), 0, 0, 0, 87.6);

    rotSpeedFull(idx)   = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
    pitchFull(idx,:)    = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);
    ptfmPitchFull(idx)  = calllib('QBladeDLL','getCustomData_at_num',ptfmPitchStr, 0, 0);
    genPwrFull(idx)     = calllib('QBladeDLL','getCustomData_at_num','Gen. Elec. Power [kW]', 0, 0);
    MytFull(idx)        = calllib('QBladeDLL','getCustomData_at_num','Y_l Mom. TWR pos 0.000 [Nm]', 0, 0);
    lssT(idx)           = calllib('QBladeDLL','getCustomData_at_num','Aero. LSS Torque [Nm]', 0, 0);
    Moop(idx)           = calllib('QBladeDLL','getCustomData_at_num','Y_c RootBend. Mom. (OOP) BLD 1 [Nm]', 0, 0);
    timeFull(idx) = calllib('QBladeDLL','getCustomData_at_num','Time [s]',0,0);    

    waitbar(idx/simTime,waitBar,'Simulation running')
end

calllib('QBladeDLL','closeInstance')
close(waitBar)

save('outputDataFinal\pi1.mat');