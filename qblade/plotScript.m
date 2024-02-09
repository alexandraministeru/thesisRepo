%% Clean environment
clearvars;clc;close all;
rng('default')

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
addpath(genpath('functions'))

%% Load data files
% Baseline
dataBaseline = load('outputDataFinal\turb_pi1.mat');

% DeePC with preview
dataDeepc = load('outputDataFinal\turb_deepc1.mat');

for idx = 1:20:length(dataDeepc.pitchFull(:,1))
    dataDeepc.pitchFull(idx,1) = dataDeepc.pitchFull(idx+1,1);
end


%% PSDs
Ts = dataBaseline.Ts;

% PSD of signals
fWaveMin = 0.05; % Hz
fWaveMax = 0.3; % Hz
fSurge = 0.0075;
fPitch = 0.0307;

% % Rotor speed
[~,~,psd_rotSpeed1] = getFFT(Ts,dataBaseline.rotSpeedFull);
[f_rotSpeed,~,psd_rotSpeed2] = getFFT(Ts,dataDeepc.rotSpeedFull);

% Platform pitch
[~,fft1,psd_ptfmPitch1] = getFFT(Ts,dataBaseline.ptfmPitchFull);
[f_ptfmPitch,fft2,psd_ptfmPitch2] = getFFT(Ts,dataDeepc.ptfmPitchFull);

% Blade pitch
[~,~,psd_bladePitch1] = getFFT(Ts,dataBaseline.pitchFull(:,1));
% [f_bladePitch,~,psd_bladePitch2] = getFFT(Ts,dataDeepc.pitchFull(:,1));
[f_bladePitch,~,psd_bladePitch2] = getFFT(1,dataDeepc.pitch_sim(:,1));

%% Plot
tsim = 0:Ts:Ts*(length(dataBaseline.rotSpeedFull)-1);
initTime = 70; %seconds
ctrlStart = 70/Ts + 1;

figure
subplot(3,2,1) % rotor speed time
plot(tsim(ctrlStart:end),dataBaseline.rotSpeedFull(ctrlStart:end),'k','LineWidth',1)
hold on
plot(tsim(ctrlStart:end),dataDeepc.rotSpeedFull(1:end-ctrlStart+1),'r','LineWidth',1)
hold on
plot(tsim(ctrlStart:end),dataDeepc.ref_y1(1:length(tsim(ctrlStart:end)))'*dataDeepc.omegaHat_max + dataDeepc.mean_rotSpeed,'k-.')
plot([tsim(ctrlStart); length(tsim)], [dataDeepc.controlParams.lby(1)*dataDeepc.omegaHat_max + dataDeepc.mean_rotSpeed; dataDeepc.controlParams.lby(1)*dataDeepc.omegaHat_max + dataDeepc.mean_rotSpeed],'k--','LineWidth',1)
plot([tsim(ctrlStart); length(tsim)], [dataDeepc.controlParams.uby(1)*dataDeepc.omegaHat_max + dataDeepc.mean_rotSpeed; dataDeepc.controlParams.uby(1)*dataDeepc.omegaHat_max + dataDeepc.mean_rotSpeed],'k--','LineWidth',1)
xlim([70 450])
% ylim([11 13.5])
ylim([8 17])
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
grid on

subplot(3,2,2) % rotor speed PSD
plot(f_rotSpeed,psd_rotSpeed1,'k','LineWidth',1)
hold on
plot(f_rotSpeed,psd_rotSpeed2,'r','LineWidth',1)
xline(fWaveMin,'k--','WaveMin','LabelVerticalAlignment','top','LineWidth',0.3)
xline(fWaveMax,'k--','WaveMax','LabelVerticalAlignment','top','LineWidth',0.3)
xline(fPitch,'k--')
xlim([0 0.4])
xlabel('Frequency (in Hz)')
ylabel('Rotor speed (in rpm^2/Hz)')
grid on

subplot(3,2,3) % platform pitch time
plot(tsim(ctrlStart:end),dataBaseline.ptfmPitchFull(ctrlStart:end),'k','LineWidth',1)
hold on
plot(tsim(ctrlStart:end),dataDeepc.ptfmPitchFull(1:end-ctrlStart+1),'r','LineWidth',1)
plot([tsim(ctrlStart); length(tsim)], [dataDeepc.controlParams.lby(2)*dataDeepc.pltfmPitchHat_max + dataDeepc.mean_ptfmPitch; dataDeepc.controlParams.lby(2)*dataDeepc.pltfmPitchHat_max + dataDeepc.mean_ptfmPitch],'k--','LineWidth',1)
plot([tsim(ctrlStart); length(tsim)], [dataDeepc.controlParams.uby(2)*dataDeepc.pltfmPitchHat_max + dataDeepc.mean_ptfmPitch; dataDeepc.controlParams.uby(2)*dataDeepc.pltfmPitchHat_max + dataDeepc.mean_ptfmPitch],'k--','LineWidth',1)
xlim([70 450])
ylim([0.5 6.5])
xlabel('Time (in s)')
ylabel('Platform pitch angle (in deg)')
grid on

subplot(3,2,4) % platform pitch PSD
plot(f_ptfmPitch,psd_ptfmPitch1,'k','LineWidth',1)
hold on
plot(f_ptfmPitch,psd_ptfmPitch2,'r','LineWidth',1)
xline(fWaveMin,'k--','WaveMin','LabelVerticalAlignment','top','LineWidth',0.3)
xline(fWaveMax,'k--','WaveMax','LabelVerticalAlignment','top','LineWidth',0.3)
xlim([0 0.4])
xlabel('Frequency (in Hz)')
ylabel('Platform pitch angle (in deg^2/Hz)')
grid on

subplot(3,2,5) % blade pitch angle time
plot(tsim(ctrlStart:end),dataBaseline.pitchFull(ctrlStart:end,1),'k','LineWidth',1)
hold on
plot(tsim(ctrlStart:end),dataDeepc.pitchFull(1:end-ctrlStart+1,1),'r','LineWidth',1)
plot([tsim(ctrlStart); length(tsim)], [dataDeepc.controlParams.lbu*dataDeepc.uhat_max + dataDeepc.mean_bladePitch; dataDeepc.controlParams.lbu*dataDeepc.uhat_max + dataDeepc.mean_bladePitch],'k--','LineWidth',1)
plot([tsim(ctrlStart); length(tsim)], [dataDeepc.controlParams.ubu*dataDeepc.uhat_max + dataDeepc.mean_bladePitch; dataDeepc.controlParams.ubu*dataDeepc.uhat_max + dataDeepc.mean_bladePitch],'k--','LineWidth',1)
xlim([70 450])
% ylim([9 14])
ylim([2.5 21.5])
xlabel('Time (in s)')
ylabel('Blade pitch angle (in deg)')
grid on

subplot(3,2,6) % blade pitch angle PSD
plot(f_bladePitch,psd_bladePitch1(1:length(f_bladePitch)),'k','LineWidth',1)
hold on
plot(f_bladePitch,psd_bladePitch2,'r','LineWidth',1)
xline(fWaveMin,'k--','WaveMin','LabelVerticalAlignment','top','LineWidth',0.3)
xline(fWaveMax,'k--','WaveMax','LabelVerticalAlignment','top','LineWidth',0.3)
xlim([0 0.4])
xlabel('Frequency (in Hz)')
ylabel('Blade pitch angle (in deg^2/Hz)')
grid on
set(gcf,'Color','White')

% %% Stats
% % LC
% LC2.rotSpeetVar(1)      = var(dataBaseline.rotSpeedFull);
% LC2.rotSpeedRMSE(1)     = rmse(dataBaseline.rotSpeedFull, 12.126*ones(length(dataBaseline.rotSpeedFull),1)); 
% LC2.bladePitchVar(1)    = var(dataBaseline.pitchFull(:,1));
% LC2.bladePitchADC(1)    = adc(dataBaseline.pitchFull(:,1),Ts,8);
% LC2.ptfmPitchVar(1)     = var(dataBaseline.ptfmPitchFull);
% genPwrFull1             = dataBaseline.rotSpeedFull*dataDeepc.ratedTq*(1/97);
% LC2.genPwrVar(1)        = var(genPwrFull1);
% LC2.MytVar(1)           = var(dataBaseline.MytFull);
% LC2.lssT(1)             = var(dataBaseline.lssT);
% LC2.Moop(1)             = var(dataBaseline.Moop);
% %%2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LC2.rotSpeetVar(2)      = var(dataDeepc.rotSpeedFull);
% LC2.rotSpeedRMSE(2)     = rmse(dataDeepc.rotSpeedFull, dataDeepc.ref_y1(1:length(dataDeepc.rotSpeedFull))'*dataDeepc.omegaHat_max + dataDeepc.mean_rotSpeed); 
% LC2.bladePitchVar(2)    = var(dataDeepc.pitchFull(:,1));
% LC2.bladePitchADC(2)    = adc(dataDeepc.pitchFull(:,1),Ts,8);
% LC2.ptfmPitchVar(2)     = var(dataDeepc.ptfmPitchFull);
% LC2.genPwrVar(2)        = var(dataDeepc.genPwrFull);
% LC2.MytVar(2)           = var(dataDeepc.MytFull);
% LC2.lssT(2)             = var(dataDeepc.lssT);
% LC2.Moop(2)             = var(dataDeepc.Moop);
% 
% % save('outputDataFinal\postProc_turb.mat','LC1')
% save('outputDataFinal\postProc_turb.mat','LC2','-append')
% 
% %% Performance
% load('outputDataFinal\postProc_turb.mat')
% 
% data = zeros(size(fieldnames(LC1),1)-3,2);
% 
% fn = fieldnames(LC1);
% for k=1:numel(fn)-3
%     data(k,1) = (LC1.(fn{k})(2)/LC1.(fn{k})(1) - 1)*100; 
% end
% 
% fn = fieldnames(LC2);
% for k=1:numel(fn)-3
%     data(k,2) = (LC2.(fn{k})(2)/LC2.(fn{k})(1) - 1)*100; 
% end
% 
% X = categorical({'Rotor speed variance','Rotor speed RMSE','Blade pitch variance','Blade pitch ADC','Platform pitch variance','Power variance'});
% X = reordercats(X,{'Power variance','Rotor speed variance','Rotor speed RMSE','Blade pitch variance','Blade pitch ADC','Platform pitch variance'});
% 
% figure
% bar(X,data)
% legend('Load case 1', 'Load case 2','Location','northwest')
% ylabel('Relative performance (%)')
% grid on
% set(gcf,'Color','White')
% 
% %% DELs
% data_del= zeros(3,1);
% Neq = 1e7;
% m  = 10;
% 
% % fn = fieldnames(LC1);
% % for k=7:numel(fn)
% %     data_del(k-6,1) = (LC1.(fn{k})(2)/LC1.(fn{k})(1) - 1)*100; 
% % end
% 
% fn = fieldnames(LC2);
% for k=7:numel(fn)
%     data_del(k-6,1) = (LC2.(fn{k})(2)/LC2.(fn{k})(1) - 1)*100; 
% end
% 
% X = categorical({'Tower base loads','LSS torque','OoP blade moment'});
% X = reordercats(X,{'Tower base loads','LSS torque','OoP blade moment'});
% 
% figure
% bar(X,data_del)
% legend('Load case 1', 'Load case 2','Location','northwest')
% ylabel('Relative DEL (%)')
% grid on
% set(gcf,'Color','White')

