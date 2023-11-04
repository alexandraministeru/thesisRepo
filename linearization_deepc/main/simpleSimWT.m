% This script runs a simulation of a WT.
%% Clear environment
clearvars;clc;close all;

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));

%% Set input FST file
OpenFASTRoot = 'D:\Master\TUD\Y2\Thesis\matlab\fromAmr\FF_OpenFAST\5MW_OC3Spar_DLL_WTurb_WavesIrr\';
FAST_InputFileName = [ OpenFASTRoot '5MW_OC3Spar_DLL_WTurb_WavesIrr.fst' ];

%% Run simulation
TMax = 3000; % seconds
sim('OpenLoop.mdl',[0,TMax]);

%% Plot nonlinear model output
outFile = 'D:\Master\TUD\Y2\Thesis\matlab\fromAmr\FF_OpenFAST\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';
plotChannels = {'RotSpeed','GenSpeed','BldPitch1','GenTq','GenPwr','Wave1Elev','B1WvsFxi','B1WvsMyi'};
PlotFASToutput(outFile,[],[],plotChannels,1)