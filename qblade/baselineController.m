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

simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_steadyWind_IrrWaves_pi.sim';
calllib('QBladeDLL','loadSimDefinition',simName)
calllib('QBladeDLL','initializeSimulation')
Ts = 0.05;
simTime = 5000;
t = 0:Ts:(simTime-1)*Ts;

rotSpeedFull2 = zeros(simTime,1);
pitchFull2 = zeros(simTime,3);
timeFull2 = zeros(simTime,3);

rotSpeedStr = 'Rotational Speed [rpm]';
pitchStr = 'Pitch Angle Blade 1 [deg]';

waitBar = waitbar(0,'Preview buffer') ;

for kBuffer = 1:1:1400
    calllib('QBladeDLL','advanceTurbineSimulation')
    calllib('QBladeDLL','advanceController_at_num', ones(1,5), 0)    
    waitbar(kBuffer/1400,waitBar,'Buffer')

end

for idx = 1:1:simTime
    calllib('QBladeDLL','advanceTurbineSimulation')
    calllib('QBladeDLL','advanceController_at_num', ones(1,5), 0)   

    rotSpeedFull2(idx) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
    pitchFull2(idx,:) = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);
    timeFull2(idx) = calllib('QBladeDLL','getCustomData_at_num','Time [s]',0,0);    

    waitbar(idx/simTime,waitBar,'Simulation running')
end

calllib('QBladeDLL','closeInstance')
close(waitBar)

save('outputData\pi.mat','rotSpeedFull2','timeFull2','pitchFull2');