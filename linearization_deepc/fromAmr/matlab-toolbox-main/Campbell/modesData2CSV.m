function filenames = modesData2CSV(BaseName, ModesData)
% Write summary of Campbell data and modes identification to several CSV files
% The files generated will be:
%    - [BaseName '_ModesID.csv ']
%    - [BaseName '_OP.csv      ']
%    - [BaseName '_PointsXX.csv'] for each operating point where XX is the index.
%
% INPUTS:
%   - BaseName    :  basename that will be used to create the different CSV Files
%   - ModesData: structure as returned by IdentifyModes
%  

nOP = length(ModesData.ModesTable);
filenames=cell(1,nOP+2);
% --- Write to individual csv files
filenames{1}=sprintf('%s_ModesID.csv'  ,BaseName);
filenames{2}=sprintf('%s_OP.csv'       ,BaseName);
cellarray2csv(filenames{1}, ModesData.modeID_table);
cellarray2csv(filenames{2}, ModesData.opTable);
for iOP =1:nOP
    filenames{iOP+2} =  sprintf('%s_Point%02d.csv',BaseName, iOP);
    cellarray2csv(filenames{iOP+2}, ModesData.ModesTable{iOP});
end
