%% Clean environment
clearvars;clc;
rng('default')

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
addpath(genpath('functions'))

%% Load data files
% Balanced DeePC
load('fullSimFinal2\steadyWind_irrWaves_Hs3Tp12_Mp_preview_QP_scaled.mat');
% load('fullSimFinal2\turbWind_irrWaves_Hs3Tp12_Fsg_Mp_preview_QP_scaled1.mat');
out1 = out;
uSeq1 = uSeq;
dataOL1 = dataOL;
uhat_max = Du(1,1);
omegaHat_max = Dy(1,1);
pltfmPitchHat_max = Dy(2,2);
pwrHat_max = Dy(3,3);

% Rotor speed DeePC
load('focusedDeePC\steadyWind_irrWaves_Hs3Tp12_Mp_preview_QP_scaled_omega.mat')
% load('focusedDeePC\turbWind_irrWaves_Hs3Tp12_Fsg_Mp_preview_QP_scaled_omega.mat')
out2 = out;
uSeq2 = uSeq;

% Platform pitch DeePC
load('focusedDeePC\steadyWind_irrWaves_Hs3Tp12_Mp_preview_QP_scaled_thetaP.mat')
% load('focusedDeePC\turbWind_irrWaves_Hs3Tp12_Fsg_Mp_preview_QP_scaled_thetaP.mat')
out3 = out;
uSeq3 = uSeq;

%% Sim data
% Get and rescale
rotSpeed1 = out1(1,:).*omegaHat_max;
ptfmPitch1 =  out1(2,:).*pltfmPitchHat_max;
bladePitch1 = rad2deg(uSeq1.*uhat_max);
genPwr1 = out1(3,:).*pwrHat_max;

rotSpeed2 = out2(1,:).*omegaHat_max;
ptfmPitch2 =  out2(2,:).*pltfmPitchHat_max;
bladePitch2 = rad2deg(uSeq2.*uhat_max);
genPwr2 = out2(3,:).*pwrHat_max;

rotSpeed3 = out3(1,:).*omegaHat_max;
ptfmPitch3 =  out3(2,:).*pltfmPitchHat_max;
bladePitch3 = rad2deg(uSeq3.*uhat_max);
genPwr3 = out3(3,:).*pwrHat_max;

%% Add back OP
% Output file of linearization
outFileOP = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.out';
plotChannels = {'RotSpeed','PtfmPitch','BldPitch1','GenPwr'};
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

rotSpeed1_full = rotSpeed1 + opVal(1);
ptfmPitch1_full = ptfmPitch1 + opVal(2);
bladePitch1_full = bladePitch1 + opVal(3);
genPwr1 = genPwr1 + opVal(4);

rotSpeed2_full = rotSpeed2 + opVal(1);
ptfmPitch2_full = ptfmPitch2 + opVal(2);
bladePitch2_full = bladePitch2 + opVal(3);
genPwr2 = genPwr2 + opVal(4);

rotSpeed3_full = rotSpeed3 + opVal(1);
ptfmPitch3_full = ptfmPitch3 + opVal(2);
bladePitch3_full = bladePitch3 + opVal(3);
genPwr3 = genPwr3 + opVal(4);

ref = 12.1*ones(length(rotSpeed1),1);

%% PSDs
% PSD of signals
fWaveMin = 0.05; % Hz
fWaveMax = 0.3; % Hz

% % Rotor speed
[~,~,psd_rotSpeed1] = getFFT(Ts,rotSpeed1);
[f_rotSpeed,~,psd_rotSpeed2] = getFFT(Ts,rotSpeed2);
[~,~,psd_rotSpeed3] = getFFT(Ts,rotSpeed3);

% Platform pitch
[~,fft1,psd_ptfmPitch1] = getFFT(Ts,ptfmPitch1);
[f_ptfmPitch,fft2,psd_ptfmPitch2] = getFFT(Ts,ptfmPitch2);
[~,fft3,psd_ptfmPitch3] = getFFT(Ts,ptfmPitch3);

% Blade pitch
[~,~,psd_bladePitch1] = getFFT(Ts,bladePitch1);
[f_bladePitch,~,psd_bladePitch2] = getFFT(Ts,bladePitch2);
[~,~,psd_bladePitch3] = getFFT(Ts,bladePitch3);

%% Plot
tsim_full = 800:1:length(rotSpeed3_full)+800-Ts;
colorVec = {'r','#26897A',[0.9290, 0.6940, 0.1250]};
figure
t = tiledlayout(3,2,'TileSpacing','tight','Padding','none');

nexttile %rotor speed time
plot(tsim_full,rotSpeed1_full,'Color',colorVec{1},'LineWidth',1)
hold on
plot(tsim_full,rotSpeed2_full,'Color',colorVec{2},'LineWidth',1)
hold on
plot(tsim_full,rotSpeed3_full,'Color',colorVec{3},'LineWidth',1)
hold on
plot(tsim_full,ref,'k--')
% xline(tsimControl(1),'k','LineWidth', 1)
% plot([tsimControl(1); length(tsim_full)], [controlParams.lby(1)*omegaHat_max + opVal(1); controlParams.lby(1)*omegaHat_max + opVal(1)],'k--','LineWidth',1)
% plot([tsimControl(1); length(tsim_full)], [controlParams.uby(1)*omegaHat_max + opVal(1); controlParams.uby(1)*omegaHat_max + opVal(1)],'k--','LineWidth',1)
xlim([800 1200])
% ylim([11 13])
% ylim([8 17])
xlabel('Time (in s)')
ylabel('Rotor speed (in rpm)')
axt1 = gca;
grid on

nexttile % rotor speed PSD
plot(f_rotSpeed,psd_rotSpeed1,'Color',colorVec{1},'LineWidth',1)
hold on
plot(f_rotSpeed,psd_rotSpeed2,'Color',colorVec{2},'LineWidth',1)
hold on
plot(f_rotSpeed,psd_rotSpeed3,'Color',colorVec{3},'LineWidth',1)
xline(fWaveMin,'k--','WaveMin','LabelVerticalAlignment','top','LineWidth',0.3)
xline(fWaveMax,'k--','WaveMax','LabelVerticalAlignment','top','LineWidth',0.3)
xlim([0 0.4])
xlabel('Frequency (in Hz)')
ylabel('Rotor speed (in rpm^2/Hz)')
axf1 = gca;
grid on

nexttile % platform pitch time
plot(tsim_full,ptfmPitch1_full,'Color',colorVec{1},'LineWidth',1)
hold on
plot(tsim_full,ptfmPitch2_full,'Color',colorVec{2},'LineWidth',1)
hold on
plot(tsim_full,ptfmPitch3_full,'Color',colorVec{3},'LineWidth',1)
% plot([tsimControl(1); length(tsim_full)], [controlParams.lby(2)*pltfmPitchHat_max + opVal(2); controlParams.lby(2)*pltfmPitchHat_max + opVal(2)],'k--','LineWidth',1)
% plot([tsimControl(1); length(tsim_full)], [controlParams.uby(2)*pltfmPitchHat_max + opVal(2); controlParams.uby(2)*pltfmPitchHat_max + opVal(2)],'k--','LineWidth',1)
% xlim('tight')
xlim([800 1200])
% ylim([1.5 4.5])
% ylim([0.5 6.5])
xlabel('Time (in s)')
ylabel('Platform pitch angle (in deg)')
axt2 = gca;
grid on

nexttile % platform pitch PSD
plot(f_ptfmPitch,psd_ptfmPitch1,'Color',colorVec{1},'LineWidth',1)
hold on
plot(f_ptfmPitch,psd_ptfmPitch2,'Color',colorVec{2},'LineWidth',1)
hold on
plot(f_ptfmPitch,psd_ptfmPitch3,'Color',colorVec{3},'LineWidth',1)
xline(fWaveMin,'k--','WaveMin','LabelVerticalAlignment','top','LineWidth',0.3)
xline(fWaveMax,'k--','WaveMax','LabelVerticalAlignment','top','LineWidth',0.3)
xlim([0 0.4])
xlabel('Frequency (in Hz)')
ylabel('Platform pitch angle (in deg^2/Hz)')
axf2 = gca;
grid on

nexttile % blade pitch angle time
plot(tsim_full,bladePitch1_full,'Color',colorVec{1},'LineWidth',1)
hold on
plot(tsim_full,bladePitch2_full,'Color',colorVec{2},'LineWidth',1)
hold on
plot(tsim_full,bladePitch3_full,'Color',colorVec{3},'LineWidth',1)
plot([tsim_full(1); length(tsim_full)], [rad2deg(controlParams.lbu*uhat_max) + opVal(3); rad2deg(controlParams.lbu*uhat_max) + opVal(3)],'k--','LineWidth',1)
plot([tsim_full(1); length(tsim_full)], [rad2deg(controlParams.ubu*uhat_max) + opVal(3); rad2deg(controlParams.ubu*uhat_max) + opVal(3)],'k--','LineWidth',1)
% xlim('tight')
xlim([800 1200])
% ylim([9 14])
% ylim([2.5 21.5])
xlabel('Time (in s)')
ylabel('Blade pitch angle (in deg)')
axt3 = gca;
grid on

nexttile % blade pitch angle PSD
plot(f_bladePitch,psd_bladePitch1,'Color',colorVec{1},'LineWidth',1)
hold on
plot(f_bladePitch,psd_bladePitch2,'Color',colorVec{2},'LineWidth',1)
hold on
plot(f_bladePitch,psd_bladePitch3,'Color',colorVec{3},'LineWidth',1)
xline(fWaveMin,'k--','WaveMin','LabelVerticalAlignment','top','LineWidth',0.3)
xline(fWaveMax,'k--','WaveMax','LabelVerticalAlignment','top','LineWidth',0.3)
xlim([0 0.4])
xlabel('Frequency (in Hz)')
ylabel('Blade pitch angle (in deg^2/Hz)')
axf3 = gca;
grid on

hold on
ph1 = plot(nan, nan, '.', 'MarkerSize', 30, 'Color', colorVec{1});
ph2 = plot(nan, nan, '.', 'MarkerSize', 30, 'Color', colorVec{2});
ph3 = plot(nan, nan, '.', 'MarkerSize', 30, 'Color', colorVec{3});
hold off
leg = legend([ph1 ph2 ph3],{'Balanced FF DeePC_{ }','FF DeePC-\Omega_{ }','FF DeePC-\Theta_P'},'NumColumns',3);
leg.Layout.Tile = 'south';
leg.ItemTokenSize = [20; 20; 20];
fontsize(leg,15,'points')

axt1.XTickLabels = {};
axt2.XTickLabels = {};
axt1.XLabel.String = '';
axt2.XLabel.String = '';

axf1.XTickLabels = {};
axf2.XTickLabels = {};
axf1.XLabel.String = '';
axf2.XLabel.String = '';
set(gcf,'Color','White')


%%
% Baseline
% Steady wind, Hs=3, Tp=12
outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs3_tp12\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
% outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs3_tp12_turbwind\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
[data, channels, units, headers] = ReadFASTtext(outFile);
plotChannels = {'RotSpeed','PtfmPitch','BldPitch1','GenPwr'};

id = find(ismember(channels,plotChannels{1}));
rotSpeedBaseline = data(:,id);
rotSpeed = rotSpeedBaseline(1:20:end);

id = find(ismember(channels,plotChannels{2}));
ptfmPitchBaseline = data(:,id);
ptfmPitch = ptfmPitchBaseline(1:20:end);

id = find(ismember(channels,plotChannels{3}));
bladePitchBaseline = data(:,id);
bladePitch = bladePitchBaseline(1:20:end);

id = find(ismember(channels,plotChannels{4}));
genPwrBaseline = data(:,id);
genPwr = genPwrBaseline(1:20:end);

% y_baseline = getLinWTResponse(outFile);
% rotSpeedBaseline = y_baseline(:,1);
% rotSpeed3 = rotSpeedBaseline(1:20:end);
% 
% ptfmPitchBaseline = y_baseline(:,2);
% ptfmPitch3 = ptfmPitchBaseline(1:20:end);
% 
% id = find(ismember(channels,plotChannels{3}));
% bladePitchBaseline = data(:,id);
% bladePitch3 = bladePitchBaseline(1:20:end);
% 
% genPwrBaseline = y_baseline(:,3);
% genPwr3 = genPwrBaseline(1:20:end);
%% Performance
data_performance = zeros(6,3);

data_performance(1,1) = (var(rotSpeed1)/var(rotSpeed) - 1)*100;
data_performance(1,2) = (var(rotSpeed2)/var(rotSpeed) - 1)*100;
data_performance(1,3) = (var(rotSpeed3)/var(rotSpeed) - 1)*100;

data_performance(2,1) = (rmse(rotSpeed1 + opVal(1),ref')/rmse(rotSpeed(800:end-2), ref) - 1)*100;
data_performance(2,2) = (rmse(rotSpeed2 + opVal(1),ref')/rmse(rotSpeed(800:end-2), ref) - 1)*100;
data_performance(2,3) = (rmse(rotSpeed3 + opVal(1),ref')/rmse(rotSpeed(800:end-2), ref) - 1)*100;

data_performance(3,1) = (var(bladePitch1(50:end))/var(bladePitch3) - 1)*100;
data_performance(3,2) = (var(bladePitch2(50:end))/var(bladePitch3) - 1)*100;
data_performance(3,3) = (var(bladePitch3(50:end))/var(bladePitch3) - 1)*100;

data_performance(4,1) = (adc(bladePitch1,Ts,8)/adc(bladePitchBaseline,0.05,8) - 1)*100;
data_performance(4,2) = (adc(bladePitch2,Ts,8)/adc(bladePitchBaseline,0.05,8) - 1)*100;
data_performance(4,3) = (adc(bladePitch3,Ts,8)/adc(bladePitchBaseline,0.05,8) - 1)*100;

data_performance(5,1) = (var(ptfmPitch1)/var(ptfmPitch) - 1)*100;
data_performance(5,2) = (var(ptfmPitch2)/var(ptfmPitch) - 1)*100;
data_performance(5,3) = (var(ptfmPitch3)/var(ptfmPitch) - 1)*100;

data_performance(6,1) = (var(genPwr1)/var(genPwr) - 1)*100;
data_performance(6,2) = (var(genPwr2)/var(genPwr) - 1)*100;
data_performance(6,3) = (var(genPwr3)/var(genPwr) - 1)*100;

%% Plot performance
X = categorical({'Rotor speed variance','Rotor speed RMSE','Blade pitch variance','Blade pitch ADC','Platform pitch variance','Power variance'});
X = reordercats(X,{'Power variance','Rotor speed variance','Rotor speed RMSE','Blade pitch variance','Blade pitch ADC','Platform pitch variance'});

figure
t = tiledlayout(1,1,'TileSpacing','tight','Padding','none');
nexttile
b= bar(X,data_performance);
b(1).FaceColor = colorVec{1};
b(2).FaceColor = colorVec{2};
b(3).FaceColor = colorVec{3};
legend('Balanced FF DeePC','FF DeePC-\Omega','FF DeePC-\Theta_P','Location','northwest')
ylabel('Relative performance (%)')
grid on
set(gcf,'Color','White')

%% DELs
data_del = zeros(3,3);
y_baseline = getLinWTResponse(outFile);
twr = y_baseline(:,4);
lss = y_baseline(:,5);
oop = y_baseline(:,6);

data_del(1,1) = (var(out1(4,:).*Dy(4,4))/var(twr) - 1)*100;
data_del(1,2) = (var(out2(4,:).*Dy(4,4))/var(twr) - 1)*100;
data_del(1,3) = (var(out3(4,:).*Dy(4,4))/var(twr) - 1)*100;

data_del(2,1) = (var(out1(5,:).*Dy(5,5))/var(lss) - 1)*100;
data_del(2,2) = (var(out2(5,:).*Dy(5,5))/var(lss) - 1)*100;
data_del(2,3) = (var(out3(5,:).*Dy(5,5))/var(lss) - 1)*100;

data_del(3,1) = (var(out1(6,:).*Dy(6,6))/var(oop) - 1)*100;
data_del(3,2) = (var(out2(6,:).*Dy(6,6))/var(oop) - 1)*100;
data_del(3,3) = (var(out3(6,:).*Dy(6,6))/var(oop) - 1)*100;

%% Plot DEL
X = categorical({'Tower base loads','LSS torque','OoP blade moment'});
X = reordercats(X,{'Tower base loads','LSS torque','OoP blade moment'});

figure
t = tiledlayout(1,1,'TileSpacing','tight','Padding','none');
nexttile
b= bar(X,data_del);
b(1).FaceColor = colorVec{1};
b(2).FaceColor = colorVec{2};
b(3).FaceColor = colorVec{3};
legend('Balanced FF DeePC','FF DeePC-\Omega','FF DeePC-\Theta_P','Location','northwest')
ylabel('Relative DEL (%)')
grid on
set(gcf,'Color','White')


figure
plot(diff(uSeq3,2))
hold on
plot(diff(uSeq1,2))
grid on


figure
plot(out2(4,:).*Dy(4,4))
hold on
plot(out3(4,:).*Dy(4,4))
grid on






