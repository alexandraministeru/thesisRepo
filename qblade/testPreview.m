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

%%
calllib('QBladeDLL','createInstance',0,24)
calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')

% projectFile = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL5MW.qpr';
% calllib('QBladeDLL','loadProject',projectFile)

simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL5MW\NREL_5MW_turbWind.sim';
calllib('QBladeDLL','loadSimDefinition',simName)

calllib('QBladeDLL','initializeSimulation')

f = 100;
avgWindSpeed = 16;

simTime = f; %in timestep, actual time is timestep*#timesteps
t = 0:Ts:(simTime-1)*Ts;

V_hub = zeros(simTime,3);
V_hub_preview = zeros(simTime,3);

for i = 1:1:simTime     
    calllib('QBladeDLL','advanceTurbineSimulation')
    if i == 1
        for idxPreview = 1:1:f
            xDist = idxPreview*Ts*avgWindSpeed;
            V_hub_preview(idxPreview,:) = calllib('QBladeDLL','getWindspeed', -xDist, 0, 87.6, [0 0 0]);
        end
    end    
    V_hub(i,:) = calllib('QBladeDLL','getWindspeed', 0, 0, 87.6, [0 0 0]);    
end
calllib('QBladeDLL','closeInstance')


figure
plot(V_hub(2:end,1))
hold on
plot(V_hub_preview(1:end-1,1))
grid on
xlabel('Time (in s)')
ylabel('Horizontal wind speed (in m/s)')
legend('Wind Speed at hub','Wind Speed preview at hub','Location','southeast')
set(gcf,'Color','White')

