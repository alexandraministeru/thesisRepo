function [LTIsys, MBC, matData, FAST_linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath)

%FASTLinearization(FilePath,MtlbTlbxPath) A funtion used to return an LTI system from OpenFAST
%
% Input arguments
%----------------
% FilePath    : Complete path to the directory containing the *.lin files.
% MtlbTlbxPath: Complete path to the directory containing OpenFAST MATLAB-Toolbox.
%
% Output arguments: 
%------------------
% LTIsys      : An LTI system in a state-space formulation.
%
%=================================================================================
% Author: Amr Hegazy
% Date  : 06-Jul-2022
% Email : <a href="mailto:a.r.hegazy@tudelft.nl">a.r.hegazy@tudelft.nl</a>
% Delft Center for Systems and Control (DCSC)
% TU Delft, The Netherlands

%=================================================================================

addpath((genpath(MtlbTlbxPath)))
%cd(FilePath)

%% Collecting the *.lin files in a cell array for the MBC transformation
FileNames = LinFileNames(FilePath,MtlbTlbxPath);

%% Running the MBC transformation routine
[MBC, matData, FAST_linData, VTK] = fx_mbc3(FileNames); 

%% Finding the location of the generator azimuth state in the state vector
GeAz_index = find(cellfun(@isempty,strfind(MBC.DescStates,'ED Variable speed generator DOF (internal DOF index = DOF_GeAz), rad')) == false);

%% Removing the generator azimuth state from both the state vector and the averaged state-space matrices
MBC.DescStates(GeAz_index) = [];

MBC.AvgA(GeAz_index,:) = [];
MBC.AvgA(:,GeAz_index) = [];
MBC.AvgB(GeAz_index,:) = [];
MBC.AvgC(:,GeAz_index) = [];

%% Finding the location of the AeroDyn UserProp in the input vector
ADup_index = find(cellfun(@isempty,strfind(MBC.DescCntrlInpt,'AD User property on blade')) == false);

%% Removing the AeroDyn UserProp inputs from both the input vector and the averaged state-space matrices
MBC.DescCntrlInpt(ADup_index) = [];
MBC.AvgB(:,ADup_index) = [];
MBC.AvgD(:,ADup_index) = [];

%% Generating an LTI state-space model
LTIsys = ss(MBC.AvgA,MBC.AvgB,MBC.AvgC,MBC.AvgD);

end
