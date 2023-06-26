function FileNames = LinFileNames(FilePath,MtlbTlbxPath)

% A routine used to collect the *.lin files in a cell array to be provided to the MBC post-processor.
%
% Input arguments
%----------------
% FilePath    : Complete path to the directory containing the *.lin files.
% MtlbTlbxPath: Complete path to the directory containing OpenFAST MATLAB-Toolbox.
%
% Output arguments
%-----------------
% FileNames: A cell array comprising all the *.lin file in the directory.
%
%=====================================================================================================
%
% Author: Amr Hegazy
% Date  : 24/July/2022
% Delft Center for Systems and Control (DCSC)
% Delft, The Netherlands

%=====================================================================================================

addpath((genpath(MtlbTlbxPath)))

cd(FilePath)

LinFiles = dir('*.lin');                    % Finding all the linearization output, (*.lin), files 
NLinTimes = length(LinFiles);               % Obtaining the number of the linearization output, (*.lin), files
FileNames = cell(NLinTimes,1);              % Pre-allocation of the cell array to be filled with the string of each linearization output, (*.lin), file

%% Filling the cell-array with the strings of the linearization output files
for i = 1:NLinTimes
    FileNames{i} = LinFiles(i).name;
end

end