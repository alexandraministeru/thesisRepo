%% Clean environment
clearvars;clc;close all;
rng('default')

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
addpath(genpath('functions'))

%% Get preview data
load('inputData\waveForces.mat');

%% % Normalized correlation of wave forces
% [c,lags] = xcorr(F_surge, M_pitch, 'coeff');
% Ts = 0.05;
% 
% figure
% stem(lags*Ts,c,'k','filled')
% hold on
% plot(0,c(lags==0),'ro','MarkerSize',9,'MarkerFaceColor','r')
% xlim([0 30])
% ylim([-1 1])
% xlabel('Time lag (in s)')
% ylabel('Correlation')
% title('Normalized correlation between F_{surge} and M_{pitch}')
% grid on
% set(gcf,'Color','White')

%% Tower loads
outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs3_tp12\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
[data, channels, units, headers] = ReadFASTtext(outFile);
plotChannels = {'TwrBsMyt'};
% id = find(ismember(channels,plotChannels{1}));
% Myt = data(:,id);
Ts = 0.05;
ssWindowIdx = 60/Ts;
opVal_Myt = getSSMean(data, ssWindowIdx, channels, plotChannels{1});

%% Load linearization
load('inputData\linDataWave.mat');

%% Find index of blade pitch, gen speed and wind speed
inputChannelsList = MBC.DescCntrlInpt;
outputChannelsList = MBC.DescOutput;

inputChannels = {'ED Extended input: collective blade-pitch command, rad', ...
    'ED Generator torque, Nm', ...
    'ED Platform X force, node 1, N', ...
    'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s',...    
    'ED Platform Y moment, node 1, Nm'};    

outputChannels = {'ED RotSpeed, (rpm)', ...
    'ED PtfmPitch, (deg)', ...
    'ED TwrBsMyt, (kN-m)'};

LTIsys_reduced = getReducedSS(MBC,LTIsys,inputChannels,outputChannels);

% Set input, output and state names
LTIsys_reduced.InputName = {'Collective blade pitch (\theta_c)', ...
    'Generator torque (\tau_g)', ...
    'Wave surge force (F_{x,P})', ...
    'Horizontal wind speed (v)', ...    
    'Wave pitch moment (M_{y,P})'}; 

LTIsys_reduced.InputUnit = {'rad','Nm','N','m/s','Nm'};
LTIsys_reduced.OutputName = {'Rotor speed (\Omega)','Platform pitch angle','Tower bending moment'};
LTIsys_reduced.OutputUnit = {'rpm','deg','kN-m'};
 
% LTIsys_reduced.StateName = {'Platform horizontal surge translation', ...
% 'Platform pitch tilt rotation', ...
% '1st tower fore-aft bending', ...
% 'First time derivative of Platform horizontal surge translation', ...
% 'First time derivative of Platform pitch tilt rotation', ...
% 'First time derivative of 1st tower fore-aft bending', ...
% 'First time derivative of Variable speed generator', ...
% 'ExctnPtfmSg1' , ...
% 'ExctnPtfmSg2' , ...
% 'ExctnPtfmSg3' , ...
% 'ExctnPtfmSg4' , ...
% 'ExctnPtfmSg5' , ...
% 'ExctnPtfmSg6' , ...
% 'ExctnPtfmSg7' , ...
% 'ExctnPtfmSg8' , ...
% 'ExctnPtfmSg9' , ...
% 'ExctnPtfmSg10', ...
% 'ExctnPtfmSg11', ...
% 'ExctnPtfmSg12', ...
% 'ExctnPtfmSg13', ...
% 'ExctnPtfmSg14', ...
% 'ExctnPtfmP1'  , ...
% 'ExctnPtfmP2'  , ...
% 'ExctnPtfmP3'  , ...
% 'ExctnPtfmP4'  , ...
% 'ExctnPtfmP5'  , ...
% 'ExctnPtfmP6'  , ...
% 'ExctnPtfmP7'  , ...
% 'ExctnPtfmP8'  , ...
% 'RdtnPtfmSg1'  , ...
% 'RdtnPtfmSg2'  , ...
% 'RdtnPtfmSg3'  , ...
% 'RdtnPtfmSg4'  , ...
% 'RdtnPtfmP1'   , ...
% 'RdtnPtfmP2'   , ...
% 'RdtnPtfmP3'   , ...
% 'RdtnPtfmP4'};
% 
% LTIsys_reduced.StateUnit = {'m', 'rad', 'm', 'm/s', 'rad/s', 'm/s', 'rad/s',...
%     '','','','','','','','','','','','','','','','','','','','','','','', ...
%     '','','','','','',''};

LTIsys_reduced.Name = 'NREL 5MW linearization around 16mps';

save('reducedSS\bp_tg_fsurg_mpitch_rotsp_ptfmp_twrmyt.mat','LTIsys_reduced');
load('reducedSS\bp_tg_fsurg_mpitch_rotsp_ptfmp_twrmyt.mat');

%% Separte transfer functions
G_hat = LTIsys_reduced(1:2,1:2);
Gd_hat = LTIsys_reduced(1:2,3:5);

%% Scale transfer functions(G_hat, Gd_hat - unscaled; G, Gd scaled)
% Scaling factors
uhat_max = 10*(pi/180); % Maximum expected input (rad)
tghat_max = 0.1*43093.55; % Maximum expected generator speed error (10% of the rated value) (Nm)
vhat_max = 0.2*16; % Maximum expected wind disturbance (m/s)
MpitchHat_max = 2*max(abs(M_pitch)); % Maximum expected wave pitch moment (Nm)
FsurgeHat_max = 2*max(abs(F_surge)); % Maximum expected wave surge force (N)
% MpitchHat_max = 2*std(M_pitch); % Maximum expected wave pitch moment (Nm)
% FsurgeHat_max = 2*std(F_surge); % Maximum expected wave surge force (N)

omegaHat_max = 0.15* 12.1; % Maximum expected rotor speed error (15% around linearization OP) (rpm)
pltfmPitchHat_max = 2; % Maximum expected platform pitch angle on top of linearization OP (deg)
mytHat_max = 0.5*opVal_Myt;

Du = diag([uhat_max tghat_max]);
Dd = diag([FsurgeHat_max vhat_max MpitchHat_max]);
% Dy = diag([omegaHat_max pltfmPitchHat_max]);
Dy = diag([omegaHat_max pltfmPitchHat_max]);% mytHat_max]);


% Scaled transfer functions
G = Dy\G_hat*Du;
Gd = Dy\Gd_hat*Dd;

% % Set input, output and state names for G
% G.InputName = {'Collective blade pitch (\theta_c)', ...
%     'Generator torque (\tau_g)'}; 
% G.InputUnit = {'-','-'};
% % G.OutputName = {'Rotor speed (\Omega)','Platform pitch angle'};
% % G.OutputUnit = {'-','-'};
% G.OutputName = {'Rotor speed (\Omega)','Platform pitch angle (\Theta_p)', 'Tower bending moment (M_{y,T})'};
% G.OutputUnit = {'-','-','-'}; 
% 
% % Set input, output and state names for Gd
% Gd.InputName = {'Wave surge force (F_{x,P})', ...
%     'Horizontal wind speed (v)', ...
%     'Wave pitch moment (M_{y,P})'}; 
% Gd.InputUnit = {'-','-','-'};
% % Gd.OutputName = {'Rotor speed (\Omega)','Platform pitch angle'};
% % Gd.OutputUnit = {'-','-'};
% Gd.OutputName = {'Rotor speed (\Omega)','Platform pitch angle (\Theta_p)', 'Tower bending moment (M_{y,T})'};
% Gd.OutputUnit = {'-','-','-'}; 

%% Bode plots
fWaveMin = 0.05;
fWaveMax = 0.3;
fSurge = 0.0075;
fPitch = 0.0307;

figure
bPlot1 = bodeplot(G);
setoptions(bPlot1,'FreqUnits','Hz','PhaseVisible','off','XLim',[1e-3,1])
h = get(gcf,'Children');
for idxSP = 3:2:13
    xline(h(idxSP),fWaveMin,'k--','WaveMin','LabelVerticalAlignment','bottom','LineWidth',0.3)
    xline(h(idxSP),fWaveMax,'k--','WaveMax','LabelVerticalAlignment','bottom','LineWidth',0.3)
    xline(h(idxSP),fSurge,'k--','Surge','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','LineWidth',0.3)
    xline(h(idxSP),fPitch,'k--','Pitch','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','LineWidth',0.3)
end
grid on
set(gcf,'Color','White')

figure
bPlot1 = bodeplot(Gd(:,1:2));
setoptions(bPlot1,'FreqUnits','Hz','PhaseVisible','off','XLim',[1e-3,1])
hold on
bPlot2 = bodeplot(Gd(:,3));
setoptions(bPlot2,'FreqUnits','Hz','PhaseVisible','off','XLim',[1e-3,1])
h = get(gcf,'Children');
for idxSP = 3:2:13
    xline(h(idxSP),fWaveMin,'k--','WaveMin','LabelVerticalAlignment','bottom','LineWidth',0.3)
    xline(h(idxSP),fWaveMax,'k--','WaveMax','LabelVerticalAlignment','bottom','LineWidth',0.3)
    xline(h(idxSP),fSurge,'k--','Surge','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','LineWidth',0.3)
    xline(h(idxSP),fPitch,'k--','Pitch','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','LineWidth',0.3)
end
grid on
set(gcf,'Color','White')

%% Conditioning number, RGA, SVD

% G_freq = evalfr(tf(G),(1/12))
% RGA1 = G_freq.*pinv(G_freq).'
% [U,singular_values,V] = svd(G_freq)
% condition_number1 = max(diag(singular_values))/min(diag(singular_values))
freqs = logspace(-2,0);
RGAvec = zeros(2,2,length(freqs));
cnVec = zeros(length(freqs),1);
sgVals = zeros(2,2,length(freqs));

for idxFreqs = 1:length(freqs)
    G_freq = evalfr(tf(G),freqs(idxFreqs));
    RGAvec(:,:,idxFreqs) = G_freq.*pinv(G_freq).';
    [U,sgVals(:,:,idxFreqs),V] = svd(G_freq);
    cnVec(idxFreqs) = max(diag(sgVals(:,:,idxFreqs)))/min(diag(sgVals(:,:,idxFreqs)));
end

figure
semilogx(freqs,squeeze(RGAvec(1,1,:)))
hold on
semilogx(freqs,squeeze(RGAvec(1,2,:)))
xline(fWaveMin,'k--','WaveMin','LabelVerticalAlignment','bottom','LineWidth',0.3)
xline(fWaveMax,'k--','WaveMax','LabelVerticalAlignment','bottom','LineWidth',0.3)
grid on
legend('From \theta_c to \Omega', 'From \tau_g to \Omega')

figure
semilogx(freqs,cnVec)

figure
semilogx(freqs,max(diag(sgVals)))
hold on
semilogx(freqs,min(diag(sgVals)))
grid on



%% PSD waves
Ts = 0.05;
[f_Fsg,fft_Fsg,psd_Fsg] = getFFT(Ts,F_surge);
[f_Mp,fft_Mp,psd_Mp] = getFFT(Ts,M_pitch);
tsim = 0:Ts:Ts*(numel(F_surge)-1);

figure
plot(f_Fsg,psd_Fsg)
xlim([0 1])
hold on
pwelch(F_surge,[],[],[],1/Ts)
hold on
periodogram(F_surge,rectwin(length(F_surge)),length(F_surge),1/Ts)
xlim([0 1])


% 
% figure
% plot(f_Mp(find(f_Mp==0):end),psd_Mp(find(f_Mp==0):end),'k-','LineWidth',0.8)
% hold on
% plot(f_Fsg(find(f_Fsg==0):end),psd_Fsg(find(f_Fsg==0):end),'b-','LineWidth',0.8)
% xlim([0 1])
% 
% 
% figure
% sp1 = subplot(2,1,1)
% plot(f_Fsg(find(f_Fsg==0):end),fft_Fsg(find(f_Fsg==0):end))
% xlabel('Frequency (in Hz)')
% ylabel('Amplitude')
% title('FFT')
% xlim([0 1])
% grid on
% sp2 = subplot(2,1,2)
% plot(f_Fsg(find(f_Fsg==0):end),psd_Fsg(find(f_Fsg==0):end))
% xlabel('Frequency (in Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% title('PSD')
% xlim([0 1])
% grid on
% linkaxes([sp1,sp2],'x');
% sgtitle('Surge force')
% set(gcf,'Color','White')

% 
% figure
% sp1 = subplot(2,1,1)
% plot(f_Mp(find(f_Mp==0):end),fft_Mp(find(f_Mp==0):end))
% xlabel('Frequency (in Hz)')
% ylabel('Amplitude')
% title('FFT')
% xlim([0 1])
% grid on
% sp2 = subplot(2,1,2)
% plot(f_Mp(find(f_Mp==0):end),psd_Mp(find(f_Mp==0):end))
% xlabel('Frequency (in Hz)')
% ylabel('Power/Frequency (dB/Hz)')
% title('PSD')
% xlim([0 1])
% grid on
% linkaxes([sp1,sp2],'x');
% sgtitle('Pitch moment')
% set(gcf,'Color','White')
