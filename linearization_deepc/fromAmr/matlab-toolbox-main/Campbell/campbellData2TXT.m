function campbellData2TXT(TxtFilename, CampbellData, nFreqOut )
% Write summary of Campbell data (prior to mode identification) to a text file
% This summary file is useful to perform "manual" mode identification
%
% INPUTS:
%   - TxtFilename : filename to write to
%   - CampbellData: cell-array of CampbellData (one per operating point), as returned for instance by getCampbellData
%
% OPTIONAL INPUTS:
%   - nFreqOut: number of frequencies to output for each operating points. Default 20. 
%
%
if ~exist('nFreqOut','var'); nFreqOut=20; end;
nModesMax=20;

nOP = length(CampbellData);


MFreq = zeros(nModesMax,nOP);
MDamp = zeros(nModesMax,nOP);
MDesc = cell( nModesMax,nOP) ;
for iOP = 1:nOP
    CD=CampbellData{iOP};
    nModesKeep=min([length(CD.Modes), nModesMax]);
    for i = 1:nModesKeep
        DescCat = ShortModeDescr(CD,i);
        %fprintf('%8.3f ; %7.4f ; %s\n',CD.NaturalFreq_Hz(i),CD.DampingRatio(i),DescCat(1:min(120,length(DescCat))));
        MFreq(i,iOP)=CD.NaturalFreq_Hz(i);
        MDamp(i,iOP)=CD.DampingRatio(i);
        MDesc{i,iOP}=DescCat;
    end
end

%% Write a summary file, like what is outputted to screen
fid=fopen(TxtFilename,'w');
for iOP =1:nOP
    try
        WS= CampbellData{iOP}.WindSpeed;
    catch
        WS= NaN;
    end
    try
        RPM = CampbellData{iOP}.RotSpeed_rpm;
    catch
        RPM = NaN;
    end
    fprintf(fid,'------------------------------------------------------------------------\n');
    fprintf(fid,'--- OP %d - WS %.1f - RPM %.2f \n',iOP, WS, RPM);
    fprintf(fid,'------------------------------------------------------------------------\n');
    for i=1:size(MFreq,1)
        fprintf(fid, '%02d ; %8.3f ; %7.4f ; %s\n',i,MFreq(i,iOP),MDamp(i,iOP),MDesc{i,iOP});
    end
end
fclose(fid);
end

function DescCat = ShortModeDescr(CD,i)
    % Returns a shor description of the mode
    Desc      = CD.Modes{i}.DescStates(CD.Modes{i}.StateHasMaxAtThisMode);
    DescCat   = ''                                                       ;
    DescCatED = ''                                                       ;
    if length(Desc)==0
        DescCat = '' ;
        DescCatED = 'NoMax -' ;
        Desc = CD.Modes{i}.DescStates(1:8);
    end
    nBD=0;
    for iD=1:length(Desc)
        s=Desc{iD};
        s=replaceModeDescription(s);
        if Desc{iD}(1:2)=='BD'
            nBD=nBD+1;
        elseif Desc{iD}(1:2)=='ED'
            DescCatED = [s ' - ' DescCatED];
        else
            DescCat = [DescCat ' - ' s];
        end
    end
    DescCat=[DescCatED, DescCat];
    if nBD>0
        DescCat = sprintf('BD%d/%d %s',nBD,sum(CD.Modes{i}.StateHasMaxAtThisMode),DescCat);
    end
end

