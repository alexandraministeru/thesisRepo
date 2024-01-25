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

%% Load linearization
load('inputData\linDataWave.mat');

%% Find index of blade pitch, gen speed and wind speed
inputChannelsList = MBC.DescCntrlInpt;
outputChannelsList = MBC.DescOutput;

inputChannels = {'HD Extended input: wave elevation at platform ref point, m'};    

outputChannels = {'HD HydroFxi, (N)', ...
     'HD HydroMyi, (N-m)', ...    
     'ED RotSpeed, (rpm)', ...
     'ED PtfmPitch, (deg)', ...
     'ED TwrBsMyt, (kN-m)'};

LTIsys_reduced = getReducedSS(MBC,LTIsys,inputChannels,outputChannels);

% Set input, output and state names
LTIsys_reduced.InputName = {'Wave elevation (m)'}; 

LTIsys_reduced.OutputName = {'HD HydroFxi','HD HydroMyi','Rotor speed (\Omega)','Platform pitch angle','Tower bending moment'};


%% Scale transfer functions(G_hat, Gd_hat - unscaled; G, Gd scaled)
% Scaling factors
nhat_max = 3; % Maximum expected input (rad)

MpitchHat_max = max(abs(M_pitch)); % Maximum expected wave pitch moment (Nm)
FsurgeHat_max = max(abs(F_surge)); % Maximum expected wave surge force (N)

omegaHat_max = 0.15* 12.1; % Maximum expected rotor speed error (15% around linearization OP) (rpm)
pltfmPitchHat_max = 3; % Maximum expected platform pitch angle on top of linearization OP (deg)
mytHat_max = 0.5*5.311596833888427e+04;

Du = nhat_max;
Dy = diag([FsurgeHat_max MpitchHat_max omegaHat_max pltfmPitchHat_max mytHat_max]);


% Scaled transfer functions
G = Dy\LTIsys_reduced*Du;

% Set input, output and state names for G
G.InputName = {'Wave elevation (m)'}; 
G.OutputName = {'HD HydroFxi','HD HydroMyi','Rotor speed (\Omega)','Platform pitch angle','Tower bending moment'};


%% Bode plots
fWaveMin = 0.05;
fWaveMax = 0.3;
fSurge = 0.0075;
fPitch = 0.0307;

figure
bPlot1 = bodeplot(LTIsys_reduced);
setoptions(bPlot1,'FreqUnits','Hz','PhaseVisible','off','XLim',[1e-3,1])
h = get(gcf,'Children');
for idxSP = 3:2:11
    xline(h(idxSP),fWaveMin,'k--','WaveMin','LabelVerticalAlignment','bottom','LineWidth',0.3)
    xline(h(idxSP),fWaveMax,'k--','WaveMax','LabelVerticalAlignment','bottom','LineWidth',0.3)
    xline(h(idxSP),fSurge,'k--','Surge','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','LineWidth',0.3)
    xline(h(idxSP),fPitch,'k--','Pitch','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','LineWidth',0.3)
end
grid on
set(gcf,'Color','White')

