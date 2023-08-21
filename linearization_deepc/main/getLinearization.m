% This script sets parameters for linearization, runs a long simulation, 
% reaches a steady state and then linearizes around that steady
% state.
%% Clear environment
clearvars;clc;close all;

%% Set matlab-toolbox path
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));

%% Modify FST file and subfiles
oldFSTName = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr.fst';
newFSTName = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.fst';

% Derived Parameters
[templateDir, baseName, ext ] = fileparts(oldFSTName); % path, file name and file extension
if strcmp(templateDir, filesep)
    templateDir = ['.' filesep];
end

%% Read and setup path of new subfiles
% Read template FST file
FP = FAST2Matlab(oldFSTName, 2); %FP are the FST parameters, specify 2 lines of header

% Path and basename for modified files
[fdir, base,~]  = fileparts(newFSTName)       ; % Basename for subfiles
fullBase        = [fdir filesep  base]         ; % Full basename for subfiles

%% Read and modify InflowFile
% New InflowWind file
fullPathIW  = [fullBase '_IW.dat']             ; % New InflowWind file
filenameIW  = [base     '_IW.dat']             ; % New InflowWind file relative to new fst file

% Read the inflow wind file
[paramIW, templateFilenameIW] = GetFASTPar_Subfile(FP, 'InflowFile', templateDir, templateDir);

% Modify parameters of inflow wind file

% % Steady wind
paramIW_mod = SetFASTPar(paramIW    ,'WindType'  ,1); % Steady wind
paramIW_mod = SetFASTPar(paramIW_mod,'HWindSpeed',16); % Set wind speed

% % Turbulent wind
% paramIW_mod = SetFASTPar(paramIW    ,'WindType'  ,3); % Turbulent wind
% paramIW_mod = SetFASTPar(paramIW_mod,'HWindSpeed',16); % Set wind speed
% paramIW_mod = SetFASTPar(paramIW_mod,'FileName_BTS','"..\5MW_Baseline\Wind\90m_16mps.bts"');

% Write the new inflow wind file
Matlab2FAST(paramIW_mod, templateFilenameIW, fullPathIW, 2); %contains 2 header lines

%% Read and modify HydroDyn
% New HydroDyn file
fullPathHydroDyn  = [fullBase '_HydroDyn.dat']             ; % New HydroDyn file
filenameHydroDyn  = [base     '_HydroDyn.dat']             ; % New HydroDyn file relative to new fst file

% % Read the HydroDyn file
% [paramHydroDyn, templateFilenameHydroDyn] = GetFASTPar_Subfile(FP, 'HydroFile', templateDir, templateDir);
%
% % Modify parameters
% paramHydroDyn_mod = SetFASTPar(paramHydroDyn, 'WaveMod', 0); % Still water
% paramHydroDyn_mod = SetFASTPar(paramHydroDyn_mod, 'WaveHs', 3); % Wave height
% paramHydroDyn_mod = SetFASTPar(paramHydroDyn_mod, 'WaveTp', 12); % Peak-spectral period of incident waves
%
% % If these two are set to 2 it's necessary to provide state-space for the
% % excitation and radiation model
% paramHydroDyn_mod = SetFASTPar(paramHydroDyn_mod  ,'Exctnmod'  ,2);
% paramHydroDyn_mod = SetFASTPar(paramHydroDyn_mod  ,'RdtnMod'  ,2);
%
% paramHydroDyn_mod = SetFASTPar(paramHydroDyn_mod  ,'RdtnTMax'  ,4000); % not used if RdtnMod<1
% paramHydroDyn_mod = SetFASTPar(paramHydroDyn_mod  ,'RdtnDT'  ,0.05); % not used if RdtnMod<1
%
% % Write the new HydroDyn file
% Matlab2FAST(paramHydroDyn_mod, templateFilenameHydroDyn, fullPathHydroDyn, 2); %contains 2 header lines

%% Read and modify ElastoDyn
reduceDOF = 1;
if reduceDOF == 1
    % New ElastoDyn file
    fullPathElastoDyn  = [fullBase '_ElastoDyn.dat']             ; % New ElastoDyn file
    filenameElastoDyn  = [base     '_ElastoDyn.dat']             ; % New ElastoDyn file relative to new fst file

    % Read the ElastoDyn file
    [paramElastoDyn, templateFilenameElastoDyn] = GetFASTPar_Subfile(FP, 'EDFile', templateDir, templateDir);

    % Modify parameters - exclude some degrees of freedom
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn    ,'FlapDOF1', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'FlapDOF2', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'EdgeDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TeetDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'DrTrDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'YawDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TwFADOF1', 'True');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TwFADOF2', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TwSSDOF1', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TwSSDOF2', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmSgDOF', 'True');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmSwDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmHvDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmRDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmPDOF', 'True');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmYDOF', 'False');
end

% Set initial conditions
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'OoPDefl'   , 0);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'IPDefl'    , 0);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'BlPitch(1)', 11.79);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'BlPitch(2)', 11.79);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'BlPitch(3)', 11.79);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'Azimuth'   , 0);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'RotSpeed'  , 12.1);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'NacYaw'    , 0);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TTDspFA'   , 0);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TTDspSS'   , 0);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmSurge' , 15.3957);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmSway'  , 0);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmHeave' , 0);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmRoll'  , 0);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmPitch' , 3.0286);
paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmYaw'   , 0);

% Write the new ElastoDyn file
Matlab2FAST(paramElastoDyn_mod, templateFilenameElastoDyn, fullPathElastoDyn, 2);

%% Read and modify AeroDyn
fullPathAeroDyn  = [fullBase '_AeroDyn.dat']             ; % New AeroDyn file
filenameAeroDyn  = [base     '_AeroDyn.dat']             ; % New AeroDyn file relative to new fst file

% Read the AeroDyn file
[paramAeroDyn, templateFilenameAeroDyn] = GetFASTPar_Subfile(FP, 'AeroFile', templateDir, templateDir);

paramAeroDyn_mod = SetFASTPar(paramAeroDyn, 'AFAeroMod', 1);
paramAeroDyn_mod = SetFASTPar(paramAeroDyn_mod, 'FrozenWake', 'False');

% Write the new AeroDyn file
Matlab2FAST(paramAeroDyn_mod, templateFilenameAeroDyn, fullPathAeroDyn, 2);

%% Read and modify ServoDyn
fullPathServoDyn  = [fullBase '_ServoDyn.dat']             ; % New ServoDyn file
filenameServoDyn  = [base     '_ServoDyn.dat']             ; % New ServoDyn file relative to new fst file

% Read the ServoDyn file
[paramServoDyn, templateFilenameServoDyn] = GetFASTPar_Subfile(FP, 'ServoFile', templateDir, templateDir);

% For linearization
paramServoDyn_mod = SetFASTPar(paramServoDyn      ,'PCMode'  ,0);
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'YCMode'  ,0);
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'VSContrl',1);
% Following params taken from definition of NREL 5MW
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'VS_Rgn2K' ,0.0255764);
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'VS_RtGnSp',1173.7);
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'VS_RtTq'  ,43093.55);

% Write the new ServoDyn file
Matlab2FAST(paramServoDyn_mod, templateFilenameServoDyn, fullPathServoDyn, 2);

%% Extract data from the FST file and modify it
% Set a given parameter in the FST file
FP_mod = SetFASTPar(FP,'Linearize','True'); % Linearization on
FP_mod = SetFASTPar(FP_mod,'CalcSteady','True'); % Linearization on
FP_mod = SetFASTPar(FP_mod,'NLinTimes',36);
FP_mod = SetFASTPar(FP_mod,'DT',0.05);
FP_mod = SetFASTPar(FP_mod,'DT_Out',0.05);
FP_mod = SetFASTPar(FP_mod,'TrimGain',0.001);
FP_mod = SetFASTPar(FP_mod,'TrimCase',3);
FP_mod = SetFASTPar(FP_mod,'LinInputs',2);

% Change the path of subfiles to point to the newly created ones
FP_mod = SetFASTPar(FP_mod,'InflowFile',['"' filenameIW '"']);
FP_mod = SetFASTPar(FP_mod,'HydroFile',['"' filenameHydroDyn '"']);
FP_mod = SetFASTPar(FP_mod,'EDFile',['"' filenameElastoDyn '"']);
FP_mod = SetFASTPar(FP_mod,'AeroFile',['"' filenameAeroDyn '"']);
FP_mod = SetFASTPar(FP_mod,'ServoFile',['"' filenameServoDyn '"']);

%% Set simulation time
TMax = 150; % seconds
FP_mod = SetFASTPar(FP_mod,'TMax',TMax);

%% Write FST file
Matlab2FAST(FP_mod, oldFSTName, newFSTName, 2); %contains 2 header lines
fclose('all');

%% Run without linearization - get a steady state to linearize around
FAST_InputFileName = newFSTName;
sim('OpenLoop.mdl',[0,TMax]);

%% Turbulent wind data
% windTime = OutData(:,1);
% windData = OutData(:,2);

%% Plot nonlinear model output
outFile = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.out';
plotChannels = {'RotSpeed','GenSpeed','BldPitch1','GenTq','GenPwr'};
PlotFASToutput(outFile,[],[],plotChannels,1)

[data, channels, ~, ~] = ReadFASTtext(outFile);

% % Plot list of selected channels
% nPlots = length(plotChannels);
% time = data(:,1);
% figure()
% for ip = 1:nPlots
%     % Find index of plot channel within data
%     id = find(ismember(channels,plotChannels{ip}));
%     subplot(nPlots, 1, ip);
%     plot(time, data(:,id))
%     xlabel('Time (s)')
%     ylabel([channels{id} ' ' units{id}])
% end

%% Simulate linearization
FilePath = "..\5MW_OC3Spar_DLL_WTurb_WavesIrr";
MtlbTlbxPath = "D:\Master\TUD\Y2\Thesis\matlab\fromAmr\matlab-toolbox-main";
[sys, MBC, matData, linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath);

% linFile = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.1.lin';
% linData = ReadFASTLinear(linFile);
% sys = ss(linData.A, linData.B, linData.C, linData.D);

% Time vector
T = data(1:end,1);

% Construct input matrix
[~,inputs] = size(sys.B);
U = zeros(inputs,length(T));

% Horizontal wind speed
idx = find(ismember(channels,'Wind1VelX'));
U(1,:) = data(1:end,idx); %#ok<FNDSB>

% Generator torque
idx = find(ismember(channels,'GenTq'));
U(8,:) = data(1:end,idx); %#ok<FNDSB>

% Collective blade pitch - uses BldPitch1 (assumes all blades received the same command)
idx = find(ismember(channels,'BldPitch1'));
U(9,:) = data(1:end,idx);

% Extract initial condition
initState = cell2mat(linData.x_op);

% Simulate system
Y = lsim(sys,U,T,initState);

%% Plot linear model output
channelsLin = linData.y_desc;
plotChannelsLin = {'ED RotSpeed, (rpm)','ED GenSpeed, (rpm)','ED BldPitch1, (deg)'};

nPlots = length(plotChannelsLin);
for i=1:nPlots
    idx = find(ismember(channelsLin,plotChannelsLin{i}));
    figure()
    plot(T,Y(:,idx))
    grid on
    title(plotChannelsLin{i})
end

%% Validation
outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr/5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.out';
plotChannels = {'RotSpeed','GenSpeed','BldPitch1'};
[data, channels, units, headers] = ReadFASTtext(outFile);
channelsLin = linData.y_desc;
plotChannelsLin = {'ED RotSpeed, (rpm)','ED GenSpeed, (rpm)','ED BldPitch1, (deg)'};

timeSamples = length(data);
timeStep = (timeSamples - 1)/TMax;
timeWindow = 60;
ssWindowIdx = 60*timeStep;

nChannels = length(channels);
opVal = zeros(nChannels,1);
for idxCh=2:nChannels
    opVal(idxCh) = getSSMean(data, ssWindowIdx, channels, channels{idxCh});
end

% Plot list of selected channels
nPlots = length(plotChannels);
time = data(:,1);
for ip = 1:nPlots
    % Find index of plot channel within data
    id = find(ismember(channels,plotChannels{ip}));
    idLin = find(ismember(channelsLin,plotChannelsLin{ip}));
    figure()
    plot(time, data(:,id))
    hold on
    plot(time, Y(:,idLin)+opVal(id)) % add back operating point value
    xlabel('Time (s)')
    ylabel([channels{id} ' ' units{id}])
    legend('Original','Linearized','Location','SouthEast')
    grid on
end

%% Load linearization

FilePath = "..\5MW_OC3Spar_DLL_WTurb_WavesIrr";
MtlbTlbxPath = "D:\Master\TUD\Y2\Thesis\matlab\repo\linearization_deepc\fromAmr\matlab-toolbox-main";
[LTIsys, MBC, matData, FAST_linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath);
cd ..\main
