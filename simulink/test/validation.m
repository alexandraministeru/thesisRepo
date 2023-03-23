%% Clear environment
clearvars;clc

%% Set matlab-toolbox path
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));

%% Load linearization
load('lin8mps.mat')

%% Modify FST file and subfiles
oldFSTName = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr.fst';
newFSTName = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_Validation8.fst';

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

% New InflowWind file
fullPathIW  = [fullBase '_IW.dat']             ; % New InflowWind file
filenameIW  = [base     '_IW.dat']             ; % New InflowWind file relative to new fst file

% New HydroDyn file
fullPathHydroDyn  = [fullBase '_HydroDyn.dat']             ; % New HydroDyn file
filenameHydroDyn  = [base     '_HydroDyn.dat']             ; % New HydroDyn file relative to new fst file

reduceDOF = 1;
if reduceDOF == 1
    % New ElastoDyn file
    fullPathElastoDyn  = [fullBase '_ElastoDyn.dat']             ; % New ElastoDyn file
    filenameElastoDyn  = [base     '_ElastoDyn.dat']             ; % New ElastoDyn file relative to new fst file
end

%% Read and modify InflowFile
% Read the inflow wind file
[paramIW, templateFilenameIW] = GetFASTPar_Subfile(FP, 'InflowFile', templateDir, templateDir);

% Modify parameters of inflow wind file
paramIW_mod = SetFASTPar(paramIW    ,'WindType'  ,3); % binary TurbSim FF
paramIW_mod = SetFASTPar(paramIW_mod,'HWindSpeed',8); % Set wind speed

% Write the new inflow wind file
Matlab2FAST(paramIW_mod, templateFilenameIW, fullPathIW, 2); %contains 2 header lines

%% Read and modify HydroDyn
% Read the HydroDyn file
[paramHydroDyn, templateFilenameHydroDyn] = GetFASTPar_Subfile(FP, 'HydroFile', templateDir, templateDir);

% Modify parameters
paramHydroDyn_mod = SetFASTPar(paramHydroDyn, 'WaveMod', 3); % Irregular waves

% Write the new HydroDyn file
Matlab2FAST(paramHydroDyn_mod, templateFilenameHydroDyn, fullPathHydroDyn, 2); %contains 2 header lines

%% Read and modify ElastoDyn
if reduceDOF == 1
    % Read the ElastoDyn file
    [paramElastoDyn, templateFilenameElastoDyn] = GetFASTPar_Subfile(FP, 'EDFile', templateDir, templateDir);
    
    % Modify parameters - exclude some degrees of freedom
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn    , 'FlapDOF2', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod, 'DrTrDOF' , 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod, 'TwFADOF2', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod, 'TwSSDOF2', 'False');
    
    % Write the new ElastoDyn file
    Matlab2FAST(paramElastoDyn_mod, templateFilenameElastoDyn, fullPathElastoDyn, 2);
end

%% Extract data from the FST file and modify it
% Set a given parameter in the FST file
FP_mod = SetFASTPar(FP,'Linearize','False'); % no linearization

% Change the path to the inflow file to point to the newly created one
FP_mod = SetFASTPar(FP_mod,'InflowFile',['"' filenameIW '"']);

% Change the path to the HydroDyn file to point to the newly created one
FP_mod = SetFASTPar(FP_mod,'HydroFile',['"' filenameHydroDyn '"']);

if reduceDOF == 1
    % Change the path to the ElastoDyn file to point to the newly created one
    FP_mod = SetFASTPar(FP_mod,'EDFile',['"' filenameElastoDyn '"']);
end

%% Set simulation time
TMax = 500; % seconds
FP_mod = SetFASTPar(FP_mod,'TMax',TMax);

%% Write FST file
Matlab2FAST(FP_mod, oldFSTName, newFSTName, 2); %contains 2 header lines

%% Run without linearization - get a steady state to linearize around
% Paste lines 68-end of the original HydroDyn file to the modified one (also starting at line 68)
FAST_InputFileName = newFSTName;
sim('OpenLoop.mdl',[0,TMax]);

%% Simulate linear model
sys = ss(linData.A, linData.B, linData.C, linData.D);

% Time vector
T = data(1:end,1); % ! data comes from nonlinear simulation

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

%% Compare
outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr/5MW_OC3Spar_DLL_WTurb_WavesIrr_Validation8.SFunc.out';
plotChannels = {'RotSpeed','GenSpeed','BldPitch1'};
plotChannelsLin = {'ED RotSpeed, (rpm)','ED GenSpeed, (rpm)','ED BldPitch1, (deg)'};
[data, channels, units, headers] = ReadFASTtext(outFile);

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
