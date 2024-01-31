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
f = 10;
simTime = 1000;
t = 0:Ts:(simTime-1)*Ts;

mPitchPreview = zeros(f*20,1);
mPitch = zeros(f*20,1);
timeDebug = zeros(f*20,1);
waitBar = waitbar(0,'Preview buffer') ;

for kBuffer = 1:1:(f*20)
    calllib('QBladeDLL','advanceTurbineSimulation')
    calllib('QBladeDLL','getCustomData_at_num','Time [s]',0,0)

    % if mod(kBuffer-1, 20) == 0
        disp('Getting preview')
        idxPreview = kBuffer;%floor(kBuffer/20) + 1;
        mPitchPreview(idxPreview) = calllib('QBladeDLL','getCustomData_at_num','Y_g Mom. Diffraction Offset [Nm]', 0, 0);
    % end

    waitbar(kBuffer/400,waitBar,'Getting preview data')

end

for idx = 1:1:f*20
    calllib('QBladeDLL','advanceTurbineSimulation')
    calllib('QBladeDLL','getCustomData_at_num','Time [s]',0,0)

    % if mod(idx-1, 20) == 0  
        disp('Reading val')
        idxWaves = idx;%floor(idx/20) + 1;
        timeDebug(idxWaves) = calllib('QBladeDLL','getCustomData_at_num','Time [s]',0,0);
        mPitch(idxWaves) = calllib('QBladeDLL','getCustomData_at_num','Y_g Mom. Diffraction [Nm]', 0, 0);
    % end

    waitbar(idx/400,waitBar,'Getting measurement data')

end

calllib('QBladeDLL','closeInstance')
close(waitBar)

figure
plot(timeDebug,mPitchPreview,LineWidth=1)
hold on
plot(timeDebug,mPitch,LineWidth=1)
grid on
xlabel('Time (in s)')
ylabel('Platform pitch moment (in Nm)')
legend('Preview','Measurement','Location','southeast')
set(gcf,'Color','White')

diff10 = mPitchPreview-mPitch;
norm(diff20)
norm(diff10)