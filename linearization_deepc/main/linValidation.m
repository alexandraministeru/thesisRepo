% This script validates a wind turbine linearization against a simulation
% of choice.
%% Clear environment
clearvars;clc;close all;
rng('default')
%% Set matlab-toolbox path
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
addpath(genpath('functions'))
%% Set data files
% Output file of linearization
outFileOP = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.out';

% Output file of simulation to compare against

% Linearization simulation
% outFile = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.out';

% Turbulent wind, still water
% outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simTurb\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';

% Steady wind, Hs=3, Tp=12
% outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs3_tp12\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';

% Turbulent wind, Hs=3, Tp=12
outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs3_tp12_turbwind\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';

% Steady wind, Hs=4.3, Tp=10
% outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs4p3_tp10\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';

% % Turbulent wind, Hs=4.3, Tp=10
% outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\hs4p3_tp10_turbwind\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';

%% Load linearization
FilePath = "..\5MW_OC3Spar_DLL_WTurb_WavesIrr";
MtlbTlbxPath = "D:\Master\TUD\Y2\Thesis\matlab\fromAmr\linearization\matlab-toolbox-main";
[LTIsys, MBC, matData, linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath);
cd ../main

%% Eliminate inputs
inputChannelsList = MBC.DescCntrlInpt;
inputChannels = {'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s', ...
    'ED Generator torque, Nm', ...
    'ED Extended input: collective blade-pitch command, rad',...
    'ED Platform Y moment, node 1, Nm'};

nIn = length(inputChannels);
nOut = length(MBC.DescOutput);
nStates = length(LTIsys.A);

B_reduced = zeros(nStates,nIn);
D_reduced = zeros(nOut,nIn);

for idx = 1:nIn
    % Find index of input channel
    id = find(ismember(inputChannelsList,inputChannels{idx}));
    B_reduced(:,idx) = MBC.AvgB(:,id);
    D_reduced(:,idx) = MBC.AvgD(:,id);
end

LTIsys_reduced = ss(MBC.AvgA,B_reduced,MBC.AvgC,D_reduced);
% LTIsys_reduced = c2d(LTIsys_reduced,0.05);

%% Check stability
isstable(LTIsys_reduced)
pzmap(LTIsys_reduced)

disp('Only negative real part eigenvalues?'); ...
    isempty(find(real(eig(LTIsys_reduced))>=0, 1))

%% Construct input matrix
% Get nonlinear model output
[data, channels, units, headers] = ReadFASTtext(outFile);

% Time vector
T = data(1:end,1);

% Input matrix
U = zeros(nIn,length(T));

inputChannels = {'Wind1VelX',... % Horizontal wind speed (m/s)
    'GenTq',... % Generator torque (kN-m)
    'BldPitch1', ... % Collective blade pitch (assumes all blades received the same input) (deg)
    'B1WvsMyi'}; % Nm

for idxCh = 1:length(inputChannels)
    U(idxCh,:) = data(1:end,ismember(channels,inputChannels{idxCh}))';
end

%% Remove linearization point value from inputs
% Get linearization point values from the original linearization file
[dataOP, channelsOP, unitsOP, ~] = ReadFASTtext(outFileOP);

% Time vector
T_OP = dataOP(1:end,1);

% Input matrix
U_OP = zeros(nIn,length(T_OP));

for idxCh = 1:length(inputChannels)
    U_OP(idxCh,:) = dataOP(1:end,ismember(channelsOP,inputChannels{idxCh}))';
end

timeSamples = length(dataOP);
timeStep = dataOP(end,1)/(timeSamples - 1);
timeWindow = 60; % seconds
ssWindowIdx = timeWindow/timeStep;

u_opVal = mean(U_OP(:,end-ssWindowIdx:end),2);
U = U - u_opVal;

%% Fix measurement units
U(2,:) = U(2,:).*1e3; % Generator torque from kN-m to N-m
U(3,:) = deg2rad(U(3,:)); % Blade pitch from deg to rad

%% Set initial condition
XINIT = zeros(length(matData.Avgxop)-1,1);
XINIT(1) = data(1,30) - matData.Avgxop(1);
XINIT(2) = deg2rad(data(1,34)) - matData.Avgxop(2);
XINIT(3) = data(1,27) - matData.Avgxop(3);
% XINIT(1:3) = matData.Avgxop(1:3);
XINIT(4:6) = matData.Avgxop(5:7); % Set the derivative initial states
% XINIT(8:end) = matData.Avgxop(9:end);

%% Simulate linearization
% x = zeros(nStates,size(data,1)+1);
x = zeros(length(XINIT),size(data,1)+1);
% x(:,1) = XINIT;
% Y = zeros(size(data,2)-1,size(data,1));
% 
% for k=1:size(data,1)
%     x(:,k+1) = LTIsys_reduced.A*x(:,k) + LTIsys_reduced.B*U(:,k);
%     Y(:,k) = LTIsys_reduced.C*x(:,k) + LTIsys_reduced.D*U(:,k);
% end
% 
% Y = Y';

[Y,tSim,xSim] = lsim(LTIsys_reduced,U,T,XINIT);
% Y(:,35) = -Y(:,35);

%% Select plotting channels
% Noninear plot channels
% plotChannels = {'RotSpeed','PtfmPitch','GenPwr','TwrBsMyt','RotTorq','RootMyc1'};
plotChannels = {'RotSpeed','PtfmPitch','TwrBsMyt'};
% ,'BldPitch1','GenSpeed',

% Linear plot channels
channelsLin = MBC.DescOutput;
plotChannelsLin = {'ED RotSpeed, (rpm)','ED PtfmPitch, (deg)','ED TwrBsMyt, (kN-m)'};
% ,'ED BldPitch1, (deg)','ED GenSpeed, (rpm)'
%% Add back linearization point
% Find steady state operating value
nChannels = length(plotChannels);
opVal = zeros(nChannels,1);
for idxCh=1:nChannels
    opVal(idxCh) = getSSMean(dataOP, ssWindowIdx, channelsOP, plotChannels{idxCh});
end

%% Plot list of selected channels
nPlots = length(plotChannels);
% fig=figure();
% for ip = 1:nPlots
%     % Find index of plot channel within data
%     id = find(ismember(channels,plotChannels{ip}));
%     idLin = find(ismember(channelsLin,plotChannelsLin{ip}));
%     subplot(nPlots,1,ip)    
%     plot(T, data(:,id),'k','LineWidth',1)
%     hold on
%     plot(T, Y(:,idLin) + opVal(ip),'r','LineWidth',1)
%     xlabel('Time (in s)','Interpreter','latex')
%     ylabel([channels{id} ' ' units{id}],'Interpreter','latex')
%     xlim([500 1200])
%     hold on    
%     rmse(data(10000:end,id), Y(10000:end,idLin) + opVal(ip))
%     vaf(data(10000:end,id), Y(10000:end,idLin) + opVal(ip)) 
% end
% ph1 = plot(nan, nan, 'k.', 'MarkerSize', 20);
% ph2 = plot(nan, nan, 'r.', 'MarkerSize', 20);
% hold off
% leg = legend([ph1 ph2],{'Nonlinear model','Linearization'},'NumColumns',2,'Interpreter','latex','Location','south');
% leg.ItemTokenSize = [10; 10];
% grid on
% set(gcf,'Color','White')

% exportgraphics(fig,'Figures\linValidation_steady.pdf','Resolution',500)

%%
figure
t = tiledlayout(nPlots,1,'TileSpacing','tight','Padding','none');
for ip = 1:nPlots
    nexttile
    id = find(ismember(channels,plotChannels{ip}));
    idLin = find(ismember(channelsLin,plotChannelsLin{ip}));
    plot(T, data(:,id),'k','LineWidth',1)
    hold on
    plot(T, Y(:,idLin) + opVal(ip),'r','LineWidth',1)
    xlabel('Time (in s)')
    ylabel([channels{id} ' ' units{id}])
    xlim([500 1000])
    grid on
    eval(['axt' num2str(ip) ' = gca;'])
end
hold on
ph1 = plot(nan, nan, 'k.', 'MarkerSize', 30);
ph2 = plot(nan, nan, 'r.', 'MarkerSize', 30);
hold off
leg = legend([ph1 ph2],{'Nonlinear model','Linearization'},'NumColumns',2);
leg.Layout.Tile = 'south';
leg.ItemTokenSize = [20; 20];
fontsize(leg,15,'points')
set(gcf,'Color','White')

% axt1.XTickLabels = {};
% axt2.XTickLabels = {};
% axt3.XTickLabels = {};
% axt1.XLabel.String = '';
% axt2.XLabel.String = '';
% axt3.XLabel.String = '';

exportgraphics(t,'D:\Master\TUD\Y2\Thesis\Report\Figures2\linearization\case3_noPadding.eps')



