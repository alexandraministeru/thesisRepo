% make sure the OpenFAST directory where the FAST_SFunc.mex* file is located
% is in the MATLAB path (also make sure any other OpenFAST library files that
% are needed are on the MATLAB path)
%    (relative path names are not recommended in addpath()):
% addpath('../../../build/bin'); % install location for Windows Visual Studio builds
% addpath(genpath('../../../install')); % cmake default install location
%% Clear environment
clearvars;clc

%% Set matlab-toolbox path
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));

%% Modify FST file and subfiles
oldFSTName = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr.fst';
newFSTName = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_Modified8.fst';

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

reduceDOF = 0;
if reduceDOF == 1
    % New ElastoDyn file
    fullPathElastoDyn  = [fullBase '_ElastoDyn.dat']             ; % New ElastoDyn file
    filenameElastoDyn  = [base     '_ElastoDyn.dat']             ; % New ElastoDyn file relative to new fst file
end

%% Read and modify InflowFile
% Read the inflow wind file
[paramIW, templateFilenameIW] = GetFASTPar_Subfile(FP, 'InflowFile', templateDir, templateDir);

% Modify parameters of inflow wind file
paramIW_mod = SetFASTPar(paramIW    ,'WindType'  ,1); % Steady wind
paramIW_mod = SetFASTPar(paramIW_mod,'HWindSpeed',8); % Set wind speed

% Write the new inflow wind file
Matlab2FAST(paramIW_mod, templateFilenameIW, fullPathIW, 2); %contains 2 header lines

%% Read and modify HydroDyn
% Read the HydroDyn file
[paramHydroDyn, templateFilenameHydroDyn] = GetFASTPar_Subfile(FP, 'HydroFile', templateDir, templateDir);

% Modify parameters
paramHydroDyn_mod = SetFASTPar(paramHydroDyn, 'WaveMod', 0); % Still water

% Write the new HydroDyn file
Matlab2FAST(paramHydroDyn_mod, templateFilenameHydroDyn, fullPathHydroDyn, 2); %contains 2 header lines

%% Read and modify ElastoDyn
if reduceDOF == 1
    % Read the ElastoDyn file
    [paramElastoDyn, templateFilenameElastoDyn] = GetFASTPar_Subfile(FP, 'EDFile', templateDir, templateDir);

    % Modify parameters - exclude some degrees of freedom
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn    ,'FlapDOF2', 'False');
    paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'DrTrDOF', 'False');
    paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'TwFADOF2', 'False');
    paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'TwSSDOF2', 'False');

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
TMax = 300; % seconds
FP_mod = SetFASTPar(FP,'TMax',TMax); % Linearization

%% Write FST file
Matlab2FAST(FP_mod, oldFSTName, newFSTName, 2); %contains 2 header lines

%% Run without linearization - get a steady state to linearize around
% Paste lines 68-end of the original HydroDyn file to the modified one (also starting at line 68)
% addpath(genpath('D:/Master/TUD/Y2/Thesis/matlab/simulink/examples')); 
FAST_InputFileName = newFSTName;
sim('OpenLoop.mdl',[0,TMax]);

%% Set initial condition for linearization in EDFile
oldFSTName = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_Modified8.fst';
newFSTName = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_Modified8_2.fst';

[fdir, base,~]  = fileparts(newFSTName)       ; % Basename for subfiles
fullBase        = [fdir filesep  base]         ; % Full basename for subfiles

FP = FAST2Matlab(oldFSTName, 2); %FP are the FST parameters, specify 2 lines of header

% Change ElastoDyn file
fullPathElastoDyn  = [fullBase '_ElastoDyn.dat']             ; % New ElastoDyn file
filenameElastoDyn  = [base     '_ElastoDyn.dat']             ; % New ElastoDyn file relative to new fst file

% Read the ElastoDyn file
[paramElastoDyn, templateFilenameElastoDyn] = GetFASTPar_Subfile(FP, 'EDFile', templateDir, templateDir);

timeSamples = length(OutData);
timeStep = (timeSamples - 1)/TMax;
timeWindow = 60;
ssWindowIdx = 60*timeStep;

% Average over timeWindow s.t transient neglected (%OR set simulation time 10^4 seconds OR set wind steps)
OoPDefl = getSSMean(OutData, ssWindowIdx, 12);
IPDefl = getSSMean(OutData, ssWindowIdx, 13);
BlPitch1 = getSSMean(OutData, ssWindowIdx, 5);
BlPitch2 = getSSMean(OutData, ssWindowIdx, 6);
BlPitch3 = getSSMean(OutData, ssWindowIdx, 7);
Azimuth = getSSMean(OutData, ssWindowIdx, 8);
RotSpeed = getSSMean(OutData, ssWindowIdx, 9);
NacYaw = getSSMean(OutData, ssWindowIdx, 11);
TTDspFA = getSSMean(OutData, ssWindowIdx, 27);
TTDspSS = getSSMean(OutData, ssWindowIdx, 28);
PtfmSurge  = getSSMean(OutData, ssWindowIdx, 30);
PtfmSway = getSSMean(OutData, ssWindowIdx, 31);
PtfmHeave = getSSMean(OutData, ssWindowIdx, 32);
PtfmRoll = getSSMean(OutData, ssWindowIdx, 33);
PtfmPitch = getSSMean(OutData, ssWindowIdx, 34);
PtfmYaw = getSSMean(OutData, ssWindowIdx, 35);

% Set initial conditions in EDFile
paramElastoDyn_mod = SetFASTPar(paramElastoDyn    ,'OoPDefl'   , OoPDefl);
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'IPDefl'    , IPDefl); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'BlPitch(1)', BlPitch1); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'BlPitch(2)', BlPitch2); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'BlPitch(3)', BlPitch3); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'Azimuth'   , Azimuth); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'RotSpeed'  , RotSpeed); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'NacYaw'    , NacYaw); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'TTDspFA'   , TTDspFA); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'TTDspSS'   , TTDspSS); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmSurge' , PtfmSurge); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmSway'  , PtfmSway); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmHeave' , PtfmHeave); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmRoll'  , PtfmRoll); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmPitch' , PtfmPitch); 
paramElastodyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmYaw'   , PtfmYaw); 

% Write the new ElastoDyn file
Matlab2FAST(paramElastoDyn_mod, templateFilenameElastoDyn, fullPathElastoDyn, 2);

%% Modify AeroDyn file for linearization
fullPathAeroDyn  = [fullBase '_AeroDyn.dat']             ; % New AeroDyn file
filenameAeroDyn  = [base     '_AeroDyn.dat']             ; % New AeroDyn file relative to new fst file

% Read the AeroDyn file
[paramAeroDyn, templateFilenameAeroDyn] = GetFASTPar_Subfile(FP, 'AeroFile', templateDir, templateDir);

paramAeroDyn_mod = SetFASTPar(paramAeroDyn, 'AFAeroMod', 1);
% Frozen wake necessary?
% Write the new AeroDyn file
Matlab2FAST(paramAeroDyn_mod, templateFilenameAeroDyn, fullPathAeroDyn, 2);

%% Modify HydroDyn file for linearization
fullPathHydroDyn  = [fullBase '_HydroDyn.dat']             ; % New HydroDyn file
filenameHydroDyn  = [base     '_HydroDyn.dat']             ; % New HydroDyn file relative to new fst file

% Read the HydroDyn file
[paramHydroDyn, templateFilenameHydroDyn] = GetFASTPar_Subfile(FP, 'HydroFile', templateDir, templateDir);

paramHydroDyn_mod = SetFASTPar(paramHydroDyn  ,'Exctnmod'  ,0);
paramHydroDyn_mod = SetFASTPar(paramHydroDyn_mod  ,'RdtnMod'  ,0);

% Write the new HydroDyn file
Matlab2FAST(paramHydroDyn_mod, templateFilenameHydroDyn, fullPathHydroDyn, 2);

%% Modify ServoDyn file for linearization
fullPathServoDyn  = [fullBase '_ServoDyn.dat']             ; % New ServoDyn file
filenameServoDyn  = [base     '_ServoDyn.dat']             ; % New ServoDyn file relative to new fst file

% Read the ServoDyn file
[paramServoDyn, templateFilenameServoDyn] = GetFASTPar_Subfile(FP, 'ServoFile', templateDir, templateDir);

paramServoDyn_mod = SetFASTPar(paramServoDyn      ,'PCMode'  ,0);
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'YCMode'  ,0);
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'VSContrl',0);

% Write the new ServoDyn file
Matlab2FAST(paramServoDyn_mod, templateFilenameServoDyn, fullPathServoDyn, 2);

%% Write parameters in .fst file
% Linearization
FP_mod = SetFASTPar(FP,'Linearize','True'); % Linearization

% Change the path of subfiles to point to the newly created ones
FP_mod = SetFASTPar(FP_mod,'EDFile',['"' filenameElastoDyn '"']);
FP_mod = SetFASTPar(FP_mod,'AeroFile',['"' filenameAeroDyn '"']);
FP_mod = SetFASTPar(FP_mod,'HydroFile',['"' filenameHydroDyn '"']);
FP_mod = SetFASTPar(FP_mod,'ServoFile',['"' filenameServoDyn '"']);

% Write FST file
Matlab2FAST(FP_mod, oldFSTName, newFSTName, 2);

%% Run simulation with linearization on
FAST_InputFileName = newFSTName;
sim('OpenLoop.mdl',[0,TMax]);

%% Plot nonlinear model output
outFile = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_Modified8_2.SFunc.out';
plotChannels = {'RotSpeed','GenSpeed','BldPitch1','GenTq','GenPwr'};
PlotFASToutput(outFile,[],[],plotChannels,1)

[data, channels, units, headers] = ReadFASTtext(outFile);
% disp('Available channels:')
% disp(channels)
% 
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
linFile = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_Modified8_2.SFunc.2.lin';
linData = ReadFASTLinear(linFile);
sys = ss(linData.A, linData.B, linData.C, linData.D);

% Time vector
T = data(1:end,1);

% Construct input matrix
[~,inputs] = size(sys.B);
U = zeros(inputs,length(T));

% Horizontal wind speed
idx = find(ismember(channels,'Wind1VelX'));
U(1,:) = data(1:end,idx);

% Generator torque
idx = find(ismember(channels,'GenTq'));
U(8,:) = data(1:end,idx);

% Collective blade pitch - uses BldPitch1 (assumes all blades received the same command)
idx = find(ismember(channels,'BldPitch1'));
U(9,:) = data(1:end,idx);

% Extract initial condition
initState = cell2mat(linData.x_op);

% Simulate system
Y = lsim(sys,U,T,initState);

%% Plot linear model output
channelsLin = linData.y_desc;
plotChannelsLin = {'ED RotSpeed, (rpm)','ED GenSpeed, (rpm)','ED BldPitch1, (deg)'}

nPlots = length(plotChannelsLin);
for i=1:nPlots
    idx = find(ismember(channelsLin,plotChannelsLin{i}));
    figure()
    plot(T,Y(:,idx))
    grid on
    title(plotChannelsLin{i})
end

%% Validation
outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr/5MW_OC3Spar_DLL_WTurb_WavesIrr_Modified8_2.SFunc.out';
plotChannels = {'RotSpeed','GenSpeed','BldPitch1'};
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
    plot(time, Y(:,idLin))
    xlabel('Time (s)')
    ylabel([channels{id} ' ' units{id}])
    legend('Original','Linearized','Location','SouthEast')
    grid on
end
