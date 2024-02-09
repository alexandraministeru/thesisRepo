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

simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_steadyWind_IrrWaves.sim';
calllib('QBladeDLL','loadSimDefinition',simName)
calllib('QBladeDLL','initializeSimulation')
Ts = 0.05;
f = 20;
% simTime = 2000;
% t = 0:Ts:(simTime-1)*Ts;
% 
% mPitchPreview = zeros(f*20,1);
% mPitch = zeros(f*20,1);
% timeDebug = zeros(simTime,1);
% waitBar = waitbar(0,'Preview buffer') ;
% 
% mPitchPreview = zeros(simTime+400,1);
% mPitch = zeros(simTime,1);
% 
% for kBuffer = 1:1:400
%     calllib('QBladeDLL','advanceTurbineSimulation')
%     calllib('QBladeDLL','getCustomData_at_num','Time [s]',0,0)
% 
%     disp('Getting preview')
%     idxPreview = kBuffer;
%     mPitchPreview(idxPreview) = calllib('QBladeDLL','getCustomData_at_num','Y_g Mom. Diffraction Offset [Nm]', 0, 0);
% 
%     waitbar(kBuffer/400,waitBar,'Getting preview data')
% 
% end
% 
% for idx = 1:1:simTime
%     calllib('QBladeDLL','advanceTurbineSimulation')
%     calllib('QBladeDLL','getCustomData_at_num','Time [s]',0,0)
% 
%     timeDebug(idx) = calllib('QBladeDLL','getCustomData_at_num','Time [s]',0,0);
%     mPitch(idx) = calllib('QBladeDLL','getCustomData_at_num','Y_g Mom. Diffraction [Nm]', 0, 0);
%     mPitchPreview(idxPreview+idx) = calllib('QBladeDLL','getCustomData_at_num','Y_g Mom. Diffraction Offset [Nm]', 0, 0);
%     waitbar(idx/simTime,waitBar,'Getting measurement data')
% 
% end
% 
% calllib('QBladeDLL','closeInstance')
% close(waitBar)
% 
% figure
% plot(timeDebug,mPitchPreview(401:end),LineWidth=1)
% hold on
% plot(timeDebug,mPitch,LineWidth=1)
% grid on
% xlabel('Time (in s)')
% ylabel('Platform pitch moment (in Nm)')
% legend('Preview','Measurement','Location','southeast')
% set(gcf,'Color','White')
% 
% rmse(mPitch,mPitchPreview(401:end))
% diff10 = mPitchPreview(401:end)-mPitch;
% norm(diff10)
% 
% timeDebugP = 0:0.05:simTime*Ts-Ts;

%%
simTime = 16000;
t = 0:Ts:(simTime-1)*Ts;
mPitch = zeros(simTime,1);
waitBar = waitbar(0,'Preview buffer') ;

for idx = 1:1:simTime
    idx
    calllib('QBladeDLL','advanceTurbineSimulation')
    mPitch(idx) = calllib('QBladeDLL','getCustomData_at_num','Y_g Mom. Diffraction [Nm]', 0, 0);
    waitbar(idx/simTime,waitBar,'Getting measurement data')

end

calllib('QBladeDLL','closeInstance')
close(waitBar)

save('inputData\waveData_long.mat')

