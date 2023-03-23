% This script sets parameters for linearization, runs a long simulation one
% time, reaches a steady state and then linearizes around that steady
% state. The linearization is simulated and the output is compared against 
% the nonlinear system.
%% Clear environment
clearvars;clc

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
paramIW_mod = SetFASTPar(paramIW    ,'WindType'  ,1); % Steady wind
paramIW_mod = SetFASTPar(paramIW_mod,'HWindSpeed',8); % Set wind speed

% Write the new inflow wind file
Matlab2FAST(paramIW_mod, templateFilenameIW, fullPathIW, 2); %contains 2 header lines

%% Read and modify HydroDyn
% New HydroDyn file
fullPathHydroDyn  = [fullBase '_HydroDyn.dat']             ; % New HydroDyn file
filenameHydroDyn  = [base     '_HydroDyn.dat']             ; % New HydroDyn file relative to new fst file

% Read the HydroDyn file
[paramHydroDyn, templateFilenameHydroDyn] = GetFASTPar_Subfile(FP, 'HydroFile', templateDir, templateDir);

% Modify parameters
paramHydroDyn_mod = SetFASTPar(paramHydroDyn, 'WaveMod', 0); % Still water
% For linearization:
paramHydroDyn_mod = SetFASTPar(paramHydroDyn_mod  ,'Exctnmod'  ,0);
paramHydroDyn_mod = SetFASTPar(paramHydroDyn_mod  ,'RdtnMod'  ,0);

% Write the new HydroDyn file
Matlab2FAST(paramHydroDyn_mod, templateFilenameHydroDyn, fullPathHydroDyn, 2); %contains 2 header lines

%% Read and modify ElastoDyn
reduceDOF = 1;
if reduceDOF == 1
    % New ElastoDyn file
    fullPathElastoDyn  = [fullBase '_ElastoDyn.dat']             ; % New ElastoDyn file
    filenameElastoDyn  = [base     '_ElastoDyn.dat']             ; % New ElastoDyn file relative to new fst file

    % Read the ElastoDyn file
    [paramElastoDyn, templateFilenameElastoDyn] = GetFASTPar_Subfile(FP, 'EDFile', templateDir, templateDir);

    % Modify parameters - exclude some degrees of freedom
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn    ,'FlapDOF2', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'DrTrDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TwFADOF2', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TwSSDOF2', 'False');

    % Write the new ElastoDyn file
    Matlab2FAST(paramElastoDyn_mod, templateFilenameElastoDyn, fullPathElastoDyn, 2);
end

%% Read and modify AeroDyn
fullPathAeroDyn  = [fullBase '_AeroDyn.dat']             ; % New AeroDyn file
filenameAeroDyn  = [base     '_AeroDyn.dat']             ; % New AeroDyn file relative to new fst file

% Read the AeroDyn file
[paramAeroDyn, templateFilenameAeroDyn] = GetFASTPar_Subfile(FP, 'AeroFile', templateDir, templateDir);

paramAeroDyn_mod = SetFASTPar(paramAeroDyn, 'AFAeroMod', 1);
% Frozen wake necessary for linearization?
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
FP_mod = SetFASTPar(FP_mod,'NLinTimes',1); 
FP_mod = SetFASTPar(FP_mod,'LinTimes',950); 
% FP_mod = SetFASTPar(FP_mod,'LinInputs',2); 

% Change the path of subfiles to point to the newly created ones
FP_mod = SetFASTPar(FP_mod,'InflowFile',['"' filenameIW '"']);
FP_mod = SetFASTPar(FP_mod,'HydroFile',['"' filenameHydroDyn '"']);
if reduceDOF == 1
    % Change the path to the ElastoDyn file to point to the newly created one
    FP_mod = SetFASTPar(FP_mod,'EDFile',['"' filenameElastoDyn '"']);
end
FP_mod = SetFASTPar(FP_mod,'AeroFile',['"' filenameAeroDyn '"']);
FP_mod = SetFASTPar(FP_mod,'ServoFile',['"' filenameServoDyn '"']);

%% Set simulation time
TMax = 1000; % seconds
FP_mod = SetFASTPar(FP_mod,'TMax',TMax);

%% Write FST file
Matlab2FAST(FP_mod, oldFSTName, newFSTName, 2); %contains 2 header lines

%% Run without linearization - get a steady state to linearize around
% Paste lines 68-end of the original HydroDyn file to the modified one (also starting at line 68)
FAST_InputFileName = newFSTName;
sim('OpenLoop.mdl',[0,TMax]);

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
linFile = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.1.lin';
linData = ReadFASTLinear(linFile);
sys = ss(linData.A, linData.B, linData.C, linData.D);

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

save('lin8mps.mat', 'linFile', 'linData', 'opVal')
