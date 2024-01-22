% This script sets parameters and then runs a simulation of a WT.
%% Clear environment
clearvars;clc;close all;

%% Set matlab-toolbox path
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));

%% Modify FST file and subfiles
oldFSTName = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr.fst';
newFSTName = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\5MW_OC3Spar_DLL_WTurb_WavesIrr.fst';

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
% paramIW_mod = SetFASTPar(paramIW    ,'WindType'  ,1); % Steady wind
% paramIW_mod = SetFASTPar(paramIW_mod,'HWindSpeed',16); % Set wind speed

% % Turbulent wind
paramIW_mod = SetFASTPar(paramIW    ,'WindType'  ,3); % Turbulent wind
paramIW_mod = SetFASTPar(paramIW_mod,'HWindSpeed',16); % Set wind speed
paramIW_mod = SetFASTPar(paramIW_mod,'FileName_BTS','"..\5MW_Baseline\Wind\90m_16mps.bts"');

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
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'GenDOF', 'True');%
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'YawDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TwFADOF1', 'True');%
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TwFADOF2', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TwSSDOF1', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'TwSSDOF2', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmSgDOF', 'True');%
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmSwDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmHvDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmRDOF', 'False');
    paramElastoDyn_mod = SetFASTPar(paramElastoDyn_mod,'PtfmPDOF', 'True');%
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

paramServoDyn_mod = SetFASTPar(paramServoDyn      ,'PCMode',  5);
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'YCMode'  ,5);
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'VSContrl',5);
% Following params taken from definition of NREL 5MW
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'VS_Rgn2K' ,0.0255764);
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'VS_RtGnSp',1173.7);
paramServoDyn_mod = SetFASTPar(paramServoDyn_mod  ,'VS_RtTq'  ,43093.55);

% Write the new ServoDyn file
Matlab2FAST(paramServoDyn_mod, templateFilenameServoDyn, fullPathServoDyn, 2);

%% Extract data from the FST file and modify it
% Set a given parameter in the FST file
FP_mod = SetFASTPar(FP,'Linearize','False'); % Linearization off
FP_mod = SetFASTPar(FP_mod,'DT',0.05);
FP_mod = SetFASTPar(FP_mod,'DT_Out',0.05);

% Change the path of subfiles to point to the newly created ones
FP_mod = SetFASTPar(FP_mod,'InflowFile',['"' filenameIW '"']);
FP_mod = SetFASTPar(FP_mod,'HydroFile',['"' filenameHydroDyn '"']);
FP_mod = SetFASTPar(FP_mod,'EDFile',['"' filenameElastoDyn '"']);
FP_mod = SetFASTPar(FP_mod,'AeroFile',['"' filenameAeroDyn '"']);
FP_mod = SetFASTPar(FP_mod,'ServoFile',['"' filenameServoDyn '"']);

%% Set simulation time
TMax = 2000; % seconds
FP_mod = SetFASTPar(FP_mod,'TMax',TMax);

%% Write FST file
Matlab2FAST(FP_mod, oldFSTName, newFSTName, 2); %contains 2 header lines
fclose('all');

%% Run simulation
FAST_InputFileName = newFSTName;
sim('OpenLoop.mdl',[0,TMax]);

%% Plot nonlinear model output
outFile = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr_simWave\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
plotChannels = {'RotSpeed','GenSpeed','BldPitch1','GenTq','GenPwr','Wave1Elev','B1WvsFxi','B1WvsMyi','PtfmPitch'};
PlotFASToutput(outFile,[],[],plotChannels,1)
% % 
% % % [data, channels, ~, ~] = ReadFASTtext(outFile);
% % % % Plot list of selected channels
% % % nPlots = length(plotChannels);
% % % time = data(:,1);
% % % figure()
% % % for ip = 1:nPlots
% % %     % Find index of plot channel within data
% % %     id = find(ismember(channels,plotChannels{ip}));
% % %     subplot(nPlots, 1, ip);
% % %     plot(time, data(:,id))
% % %     xlabel('Time (s)')
% % %     ylabel([channels{id} ' ' units{id}])
% % % end
