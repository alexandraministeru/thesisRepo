%% Clean environment
clearvars;clc;close all;
rng('default')

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
addpath(genpath('functions'))

%% Load data files
% DeePC no preview
load('fullSimFinal\steadyWind_irrWaves_Hs3Tp12_noPreview_QP_scaled.mat');
% load('fullSimFinal\turbWind_irrWaves_Hs3Tp12_noPreview_QP_scaled.mat')
out1 = out;
uSeq1 = uSeq;
dataOL1 = dataOL;
uhat_max = Du(1,1);
omegaHat_max = Dy(1,1);
pltfmPitchHat_max = Dy(2,2);

% DeePC with preview
load('fullSimFinal\steadyWind_irrWaves_Hs3Tp12_Mp_preview_QP_scaled.mat');
% load('fullSimFinal\turbWind_irrWaves_Hs3Tp12_Fsg_Mp_preview_QP_scaled.mat')
out2 = out;
uSeq2 = uSeq;
dataOL2 = dataOL;

% Baseline
% Steady wind, Hs=3, Tp=12
outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs3_tp12\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
% outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs3_tp12_turbwind\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
[data, channels, units, headers] = ReadFASTtext(outFile);
plotChannels = {'RotSpeed','PtfmPitch','BldPitch1'};

id = find(ismember(channels,plotChannels{1}));
rotSpeedBaseline = data(:,id);
rotSpeed3 = rotSpeedBaseline(1:20:end);

id = find(ismember(channels,plotChannels{2}));
ptfmPitchBaseline = data(:,id);
ptfmPitch3 = ptfmPitchBaseline(1:20:end);

id = find(ismember(channels,plotChannels{3}));
bladePitchBaseline = data(:,id);
bladePitch3 = bladePitchBaseline(1:20:end);


%% OL data
% Get data and rescale
tsimOL1 = dataOL1.tsim;
uOL1 = rad2deg(dataOL1.u.*uhat_max);
y1OL1 = dataOL1.y(:,1).*omegaHat_max;
y2OL1 = dataOL1.y(:,2).*pltfmPitchHat_max;

tsimOL2 = dataOL2.tsim;
uOL2 = rad2deg(dataOL2.u.*uhat_max);
y1OL2 = dataOL2.y(:,1).*omegaHat_max;
y2OL2 = dataOL2.y(:,2).*pltfmPitchHat_max;


%% Sim data
% Get and rescale
rotSpeed1 = out1(1,:).*omegaHat_max;
ptfmPitch1 =  out1(2,:).*pltfmPitchHat_max;
bladePitch1 = rad2deg(uSeq1.*uhat_max);

rotSpeed2 = out2(1,:).*omegaHat_max;
ptfmPitch2 =  out2(2,:).*pltfmPitchHat_max;
bladePitch2 = rad2deg(uSeq2.*uhat_max);


%% Concatenate
rotSpeed1_full = [y1OL1; rotSpeed1'];
ptfmPitch1_full = [y2OL1; ptfmPitch1'];
bladePitch1_full = [uOL1; bladePitch1'];

rotSpeed2_full = [y1OL2; rotSpeed2'];
ptfmPitch2_full = [y2OL2; ptfmPitch2'];
bladePitch2_full = [uOL2; bladePitch2'];

tsim_full = 0:Ts:length(rotSpeed1_full)-Ts;

%% Add back OP
% Output file of linearization
outFileOP = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.out';

[dataOP, channelsOP, unitsOP, ~] = ReadFASTtext(outFileOP);

timeSamples = length(dataOP);
timeStep = dataOP(end,1)/(timeSamples - 1);
timeWindow = 60; % seconds
ssWindowIdx = timeWindow/timeStep;

nChannels = length(plotChannels);
opVal = zeros(nChannels,1);
for idxCh=1:nChannels
    opVal(idxCh) = getSSMean(dataOP, ssWindowIdx, channelsOP, plotChannels{idxCh});
end

rotSpeed1_full = rotSpeed1_full + opVal(1);
ptfmPitch1_full = ptfmPitch1_full + opVal(2);
bladePitch1_full = bladePitch1_full + opVal(3);

rotSpeed2_full = rotSpeed2_full + opVal(1);
ptfmPitch2_full = ptfmPitch2_full + opVal(2);
bladePitch2_full = bladePitch2_full + opVal(3);

ref = 12.1*ones(length(rotSpeed1));
tsimControl = tsim + length(uOL1);

%% PSDs

% PSD of signals
fWaveMin = 0.05; % Hz
fWaveMax = 0.3; % Hz

% % Rotor speed
[~,~,psd_rotSpeed1] = getFFT(Ts,rotSpeed1);
[f_rotSpeed,~,psd_rotSpeed2] = getFFT(Ts,rotSpeed2);
[~,~,psd_rotSpeed3] = getFFT(Ts,rotSpeed3);

% Platform pitch
[~,~,psd_ptfmPitch1] = getFFT(Ts,ptfmPitch1);
[f_ptfmPitch,~,psd_ptfmPitch2] = getFFT(Ts,ptfmPitch2);
[~,~,psd_ptfmPitch3] = getFFT(Ts,ptfmPitch3);

% Blade pitch
[~,~,psd_bladePitch1] = getFFT(Ts,bladePitch1);
[f_bladePitch,~,psd_bladePitch2] = getFFT(Ts,bladePitch2);
[~,~,psd_bladePitch3] = getFFT(Ts,bladePitch3);

%% Plot
figure
subplot(3,2,1) % rotor speed time
plot(tsim_full(600:end),rotSpeed3(600:length(tsim_full)),'k','LineWidth',1)
hold on
plot(tsim_full(600:end),rotSpeed1_full(600:end),'b','LineWidth',1)
hold on
plot(tsim_full(600:end),rotSpeed2_full(600:end),'r','LineWidth',1)
hold on
plot(tsimControl,ref,'k--')
xline(tsimControl(1),'k','LineWidth', 1)
yline(controlParams.lby(1)*omegaHat_max + opVal(1),'k--','LineWidth',1)
yline(controlParams.uby(1)*omegaHat_max + opVal(1),'k--','LineWidth',1)
xlim([700 1200])
ylim([11 13])
% ylim([8 17])
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
grid on

subplot(3,2,2) % rotor speed PSD
plot(f_rotSpeed,psd_rotSpeed3(1:length(f_rotSpeed)),'k','LineWidth',1)
hold on
plot(f_rotSpeed,psd_rotSpeed1,'b','LineWidth',1)
hold on
plot(f_rotSpeed,psd_rotSpeed2,'r','LineWidth',1)
xline(fWaveMin,'k--','WaveMin','LabelVerticalAlignment','top','LineWidth',0.3)
xline(fWaveMax,'k--','WaveMax','LabelVerticalAlignment','top','LineWidth',0.3)
xlim([0 0.4])
xlabel('Frequency (in Hz)')
ylabel('Rotor speed (in rpm^2/Hz)')
grid on

subplot(3,2,3) % platform pitch time
plot(tsim_full(600:end),ptfmPitch3(600:length(tsim_full)),'k','LineWidth',1)
hold on
plot(tsim_full(600:end),ptfmPitch1_full(600:end),'b','LineWidth',1)
hold on
plot(tsim_full(600:end),ptfmPitch2_full(600:end),'r','LineWidth',1)
xline(tsimControl(1),'k','LineWidth', 1)
plot([tsimControl(1); length(tsim_full)], [controlParams.lby(2)*pltfmPitchHat_max + opVal(2); controlParams.lby(2)*pltfmPitchHat_max + opVal(2)],'k--','LineWidth',1)
plot([tsimControl(1); length(tsim_full)], [controlParams.uby(2)*pltfmPitchHat_max + opVal(2); controlParams.uby(2)*pltfmPitchHat_max + opVal(2)],'k--','LineWidth',1)
% xlim('tight')
xlim([700 1200])
ylim([1.5 4.5])
xlabel('Time (in s)')
ylabel('Platform pitch angle (in deg)')
grid on

subplot(3,2,4) % platform pitch PSD
plot(f_ptfmPitch,psd_ptfmPitch3(1:length(f_ptfmPitch)),'k','LineWidth',1)
hold on
plot(f_ptfmPitch,psd_ptfmPitch1,'b','LineWidth',1)
hold on
plot(f_ptfmPitch,psd_ptfmPitch2,'r','LineWidth',1)
xline(fWaveMin,'k--','WaveMin','LabelVerticalAlignment','top','LineWidth',0.3)
xline(fWaveMax,'k--','WaveMax','LabelVerticalAlignment','top','LineWidth',0.3)
xlim([0 0.4])
xlabel('Frequency (in Hz)')
ylabel('Platform pitch angle (in deg^2/Hz)')
grid on

subplot(3,2,5) % blade pitch angle time
plot(tsim_full(600:end),bladePitch3(600:length(tsim_full)),'k','LineWidth',1)
hold on
plot(tsim_full(600:end),bladePitch1_full(600:end),'b','LineWidth',1)
hold on
plot(tsim_full(600:end),bladePitch2_full(600:end),'r','LineWidth',1)
xline(tsimControl(1),'k','LineWidth', 1)
plot([tsimControl(1); length(tsim_full)], [rad2deg(controlParams.lbu*uhat_max) + opVal(3); rad2deg(controlParams.lbu*uhat_max) + opVal(3)],'k--','LineWidth',1)
plot([tsimControl(1); length(tsim_full)], [rad2deg(controlParams.ubu*uhat_max) + opVal(3); rad2deg(controlParams.ubu*uhat_max) + opVal(3)],'k--','LineWidth',1)
% xlim('tight')
xlim([700 1200])
ylim([9 14])
% ylim([2.5 22.5])
xlabel('Time (in s)')
ylabel('Blade pitch angle (in deg)')
grid on

subplot(3,2,6) % blade pitch angle PSD
plot(f_bladePitch,psd_bladePitch3(1:length(f_rotSpeed)),'k','LineWidth',1)
hold on
plot(f_bladePitch,psd_bladePitch1,'b','LineWidth',1)
hold on
plot(f_bladePitch,psd_bladePitch2,'r','LineWidth',1)
xline(fWaveMin,'k--','WaveMin','LabelVerticalAlignment','top','LineWidth',0.3)
xline(fWaveMax,'k--','WaveMax','LabelVerticalAlignment','top','LineWidth',0.3)
xlim([0 0.4])
xlabel('Frequency (in Hz)')
ylabel('Blade pitch angle (in deg^2/Hz)')
grid on
set(gcf,'Color','White')






