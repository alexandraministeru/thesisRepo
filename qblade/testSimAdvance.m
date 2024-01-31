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

%% Control loop
rotSpeedStr = 'Rotational Speed [rpm]';
mPitchStr = 'Y_g Mom. Diffraction [Nm]';
pitchStr = 'Pitch Angle Blade 1 [deg]';

Ts = 0.05;
tFinal = 200;
tSim = 0:Ts:tFinal-Ts;

pitch_sim = zeros(length(tSim),3);
ratedTq = 43093.55;

pitCom = idinput(length(tSim),'PRBS',[0 1/10],[-2 2]); % in degrees
pitCom = pitCom + 11.6679;
rotSpeed_out = zeros(length(tSim),1);
timeDebug = zeros(length(tSim),1);
genTorqueQB = zeros(length(tSim),1);

calllib('QBladeDLL','createInstance',0,24)
calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')

% Offshore, steady wind, irregular waves
% simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_steadyWind_IrrWaves.sim';
simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL_5MW_OC3_SPAR_steadyWind16mps_stillWater.sim';

calllib('QBladeDLL','loadSimDefinition',simName)
calllib('QBladeDLL','initializeSimulation')

waitBar = waitbar(0,'Checking advance');

for idx = 1:1:length(tSim)    
    % Advance simulation
    calllib('QBladeDLL','advanceTurbineSimulation')

    % Read inputs/outputs
    timeDebug(idx) = calllib('QBladeDLL','getCustomData_at_num','Time [s]',0,0);
    rotSpeed_out(idx) = calllib('QBladeDLL','getCustomData_at_num',rotSpeedStr, 0, 0);
    pitch_sim(idx,:) = calllib('QBladeDLL','getCustomData_at_num',pitchStr, 0, 0);
    genTorqueQB(idx) = calllib('QBladeDLL','getCustomData_at_num','Gen. HSS Torque [Nm]', 0, 0) ;

    % Set control
    if idx < 2500
        theta_c = 11;
    else
        theta_c = 11;
        ratedTq = 44093.55;
    end

    % theta_c = pitCom(idx);
    calllib('QBladeDLL','setControlVars_at_num',[ratedTq 0 theta_c theta_c theta_c],0);      
    
    waitbar(idx/length(tSim),waitBar,'checking')
end
close(waitBar)
calllib('QBladeDLL','closeInstance')
%%
% genCom = 43093.55*ones(length(tSim),1);
% genCom(2500:end) = 44093.55;
% 
pitCom = 11*ones(length(tSim),1);
pitCom(2500:end) = 15;

f = figure;
ax(1) = subplot(2,1,1);
plot(tSim,pitCom,'k','LineWidth',1)
% plot(tSim,genCom,'k','LineWidth',1)
hold on
plot(tSim,pitch_sim(:,1),'r','LineWidth',1)
% plot(tSim,genTorqueQB,'r','LineWidth',1)
xline((2500-1)*Ts,'r--')
legend('pitchCom','pitch sim')
grid on
ax(2) = subplot(2,1,2);
plot(tSim,rotSpeed_out,'k','LineWidth',1)
xline((2500-1)*Ts,'r--')
grid on
% colororder(f,"meadow")
linkaxes(ax,'x')



% %%
% calllib('QBladeDLL','createInstance',0,24)
% calllib('QBladeDLL','setLibraryPath','D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\')
% 
% % Offshore, steady wind, irregular waves
% % simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL5MW_OC3_steadyWind_IrrWaves.sim';
% simName = 'D:\Program Files\QBladeCE_2.0.6.4_win\QBladeCE_2.0.6.4\SampleProjects\NREL_5MW_OC3_SPAR_RWT\NREL_5MW_OC3_SPAR_steadyWind16mps_stillWater.sim';
% 
% calllib('QBladeDLL','loadSimDefinition',simName)
% calllib('QBladeDLL','initializeSimulation')
% simTime = 400; %in timestep, actual time is timestep*#timesteps
% valuestr = 'Rotational Speed [rpm]';
% valuestr2 = 'Gen. HSS Torque [Nm]';
% % valuestr2 = 'Gen. Power (w.o. losses) [kW]';
% 
% K = 2.24;
% N = 97;
% for i = 1:1:simTime
% 
%     calllib('QBladeDLL','advanceTurbineSimulation')
%     omega = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0.5, 0) 
%     genTorqueQB = calllib('QBladeDLL','getCustomData_at_num',valuestr2, 0, 0) 
%     genTorqueQB_store(i,:) = genTorqueQB;
% 
%     omega_g = omega*N;
%     genTorque = K.*(omega_g*(2*pi/60))^2;
%     genTorque_store(i,:) = genTorque ;
% 
%     calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 0 0 0],0)
% 
% end
% 
% calllib('QBladeDLL','closeInstance')
% 
% figure;
% plot(genTorqueQB_store)
% hold on
% plot(genTorque_store)
% grid on
% legend('QB HSS Torque','K omega^2')
% 
% 
% 


