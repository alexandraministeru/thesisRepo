function [outData]=PlotFASToutput(FASTfiles,FASTfilesDesc,ReferenceFile,Channels,ShowLegend,CustomHdr,PlotPSDs,OnePlot)
%..........................................................................
%function [timeSeriesData] = PlotFASToutput(FASTfiles,FASTfilesDesc,ReferenceFile,Channels)
%function [timeSeriesData] = PlotFASToutput(FASTfiles,FASTfilesDesc,ReferenceFile,Channels,ShowLegend,CustomHdr,PlotPSDs,OnePlot)
%
% (c) 2014-2016 National Renewable Energy Laboratory
%  Author: B. Jonkman, NREL/NWTC
%
% This routine produces time-series plots of FAST output
%..........................................................................
% Required Inputs:
% FASTfiles     - a cell array of strings, listing FAST .out or .outb file 
%                 names, whose channels are to be plotted
% Optional Inputs:
% FASTfilesDesc - a cell array of strings describing the FAST files 
%                 listed in FASTfiles, used in the plot legend. If omitted,  
%                 the routine will list them as File 1, File 2, etc.
% ReferenceFile - scalar (index into FASTfiles) that denotes which file is 
%                 considered the reference. The channels from this file 
%                 will be plotted, and the channel names from this file 
%                 will be used. If omitted, ReferenceFile is the last file 
%                 in FASTfiles. 
% Channels      - an array of channel numbers or cell array of names, 
%                 indicating which channels from the ReferenceFile
%                 should be plotted. If omitted, or is the string 'all', 
%                 all channels except time (the first one) are plotted.
% ShowLegend    - scalar logical that determines if a legend is printed.                 
% CustomHdr     - cell array describing text file format. Default will use
%                 values appropriate for FAST text output files. 
%     CustomHdr{1} = delim: delimiter for channel columns; if
%                    empty ([]), columns are delimited by whitespace
%     CustomHdr{2} = HeaderRows: number of rows in the file header, before
%                    data is encountered
%     CustomHdr{3} = NameLine: scalar value denoting line containing
%                    channel names
%     CustomHdr{4} = UnitsLine: scalar value denoting line containing
%                    channel units
% PlotPSDs      - scalar logical that determines if PSD plots will be
%                 generated in addition to the time series plots. Default
%                 is false (no PSD plots)
% OnePlot       - scalar logical that determines if each time series plot
%                 will be placed on the same or separate plots. Default
%                 is false (many plots).
%
% Note: the channels in the files need not be in the same order, but the
%  channel names must be the same [it does search for negatives]. 
%..........................................................................

if ~iscell(FASTfiles)
    FASTfiles = {FASTfiles};
end
    
numFiles = length(FASTfiles);
if numFiles < 1 
    disp('PlotFASToutput:No files to plot.')
    return
end

if nargin < 8 || isempty(OnePlot)
    OnePlot = false;
end

if nargin < 7 || isempty(PlotPSDs)
    PlotPSDs = false;
end
if nargin < 6 || isempty(CustomHdr)
    useCustomHdr=false;
else
    useCustomHdr=true;
end
if nargin < 5 || isempty(ShowLegend) 
    ShowLegend = true;
end
if nargin < 4 || isempty(Channels) 
    Channels = 'all';
end
if nargin < 3 || isempty(ReferenceFile) || (ReferenceFile < 1) || (ReferenceFile > numFiles)
    ReferenceFile = 1;
end
if nargin < 2 
    FASTfilesDesc = ''; %empty string
end


%% -----------------------------------------------------------
% Read the FAST file(s):
% ------------------------------------------------------------
data         = cell(numFiles,1);
columnTitles = cell(numFiles,1);
columnUnits  = cell(numFiles,1);
DescStr      = cell(numFiles,1);

outData      = cell(numFiles,2);

for iFile=1:numFiles

    if length(FASTfiles{iFile}) > 4 && strcmpi( FASTfiles{iFile}((end-4):end),'.outb' )
        [data{iFile}, columnTitles{iFile}, columnUnits{iFile}, ~, DescStr{iFile}] = ReadFASTbinary(FASTfiles{iFile});
    elseif ~useCustomHdr
        [data{iFile}, columnTitles{iFile}, columnUnits{iFile},    DescStr{iFile}] = ReadFASTtext(FASTfiles{iFile});                        
    else % allow other files 
        delim     = CustomHdr{1};
        HeaderRows= CustomHdr{2};
        NameLine  = CustomHdr{3};
        UnitsLine = CustomHdr{4};
        DescStr{iFile} = '';
        
        [data{iFile}, columnTitles{iFile}, columnUnits{iFile} ] = ReadFASTtext(FASTfiles{iFile}, delim, HeaderRows, NameLine, UnitsLine );    
    end
    
end

%% -----------------------------------------------------------
% Set some default values, if it wasn't input to this routine:
% ------------------------------------------------------------
if ~iscell(Channels) && strcmpi(Channels,'all')
    Channels = 2:size(data{ReferenceFile},2); % use the number of columns        
elseif ~isnumeric(Channels)    
    if iscell(Channels)
        ChannelIndx = zeros(1,length(Channels));
        for i=1:length(Channels)
            tmp = getColIndx( Channels{i}, columnTitles{ReferenceFile}, FASTfiles{ReferenceFile} );            
            if ~isempty(tmp)
                ChannelIndx(i) = tmp;
            else
                disp(['Error: ' Channels{i} ' not found in reference file.']);
            end
            
        end
        Channels = ChannelIndx(ChannelIndx ~= 0);
    else
        Channels = getColIndx( Channels, columnTitles{ReferenceFile}, FASTfiles{ReferenceFile} );
    end
end
    
if isempty(FASTfilesDesc)
    for i=1:numFiles
        FASTfilesDesc{i} = sprintf('File %i',i);
    end
end

%% -----------------------------------------------------------
% Set some default values for the plots:
% ------------------------------------------------------------
LineWidthConst = 3;
FntSz          = 17;

LineColors     = {[0 0 0],[0 1 1],[1 0 1],[0 1 0],[0 0 1],[1 0 0]};
Markers        = {'o'    ,'s'    ,'d'    ,'v'    ,'^'    ,'.'    };
% LineColors{iFile}*0.75

if numFiles > length(LineColors)
    tmp=jet(numFiles);
    for i=1:numFiles
        LineColors{i}=tmp(i,:);
    end
    n=length(Markers);
    for i=n+1:numFiles
        Markers{i}=Markers{n};
    end
end
    
titleText = char(DescStr{ReferenceFile});
FASTv8Text='Description from the FAST input file:';
indx=strfind(titleText,FASTv8Text);
if ~isempty(indx)
    titleText=titleText((indx+length(FASTv8Text)):end);
end

%% -----------------------------------------------------------
% save data for output
TimeIndex = 1;
% TimeIndex = 2;
for iFile=1:numFiles
    outData{iFile,1} = data{iFile}(:,TimeIndex);
    outData{iFile,2} = zeros( length(outData{iFile,1}), length(Channels));
end

i=0;
for iChannel = Channels
    i = i+1;    
    for iFile=1:numFiles
        
        if iChannel<1 
            err = true;
            ChannelName = ['Input Channel ' num2str(i)];
        else        
            ChannelName = columnTitles{ReferenceFile}{iChannel};       
            [ChannelIndx, err, ChannelName, scaleFact] = getColIndx( ChannelName, columnTitles{iFile}, FASTfiles{iFile} );
            
            if ~err && ChannelIndx > 0
                [columnUnits{iFile}{ChannelIndx}, scalefact2] = ConvertUnits( columnUnits{iFile}{ChannelIndx} );            
                scaleFact = scalefact2*scaleFact;
            end
        end
        
        
%         if iFile == 1
%             [ChannelName,scaleFact] = getOldFASTv8ChannelName(ChannelName);
%         elseif iFile == 2
%             [ChannelName,scaleFact] = getFASTv7ChannelName(ChannelName);
%         if iFile == numFiles
%              [ChannelName,scaleFact] = getFASTv7ChannelName(ChannelName);
%         end
        
        if err 
            outData{iFile,2}(:,i) = NaN;
            outData{iFile,3}{i}   = [ ChannelName '(not found)'];
        else
            outData{iFile,2}(:,i) = data{iFile}(:,ChannelIndx)*scaleFact;
            if scaleFact ~= 1
                outData{iFile,3}{i} = [ ChannelName ' x ' num2str(scaleFact) ];
            else
                outData{iFile,3}{i}   = ChannelName;
            end
        end
    end
end


Channels = max(Channels,1); %if any channel was missing, we'll say it's channel 1

%% -----------------------------------------------------------
% Plot the time series from each file, with each channel in 
%      a separate figure:
% ------------------------------------------------------------
[~,outFileRoot] = fileparts( FASTfiles{ReferenceFile} ); 
if OnePlot
    figNo=figure;
    ShowThisLegend = false;
else
    figNo = -1';
    ShowThisLegend = ShowLegend;
end

plotTimeSeriesData( outData, FASTfilesDesc, Markers, LineColors, ...
                    columnTitles{ReferenceFile}([TimeIndex Channels]), ...
                    columnUnits{ReferenceFile}([TimeIndex Channels]), titleText, ...
                    ShowThisLegend, LineWidthConst, FntSz, figNo, outFileRoot );
if OnePlot
    if ShowLegend
       legend show
    end
end
            

%% -----------------------------------------------------------
% Plot the psd from each file, with each channel in 
%      a separate figure:
% ------------------------------------------------------------
if PlotPSDs
    plotPSDData( outData, FASTfilesDesc, Markers, LineColors, ...
                        columnTitles{ReferenceFile}([TimeIndex Channels]), ...
                        columnUnits{ReferenceFile}([TimeIndex Channels]), titleText, ...
                        ShowLegend, LineWidthConst, FntSz );                
end

return
end
      
%%
function [f] = plotTimeSeriesData( outData, FASTfilesDesc, Markers, LineColors, ...
                RefColumnTitles, RefColumnUnits, titleText, ShowLegend, LineWidthConst, FntSz, figNo, outFileRoot )

numCols  = size(outData{1,2},2) ;
numFiles = size(outData,1);

% RefColumnTitles= columnTitles{ReferenceFile}(Channels);
% RefColumnUnits = columnUnits{ReferenceFile}(Channels);
%% -----------------------------------------------------------
% Plot the time series from each file, with each channel in 
%      a separate figure:
% ------------------------------------------------------------

    for i = 1:numCols    
        if figNo < 0
            f=figure;
            lStyle = '-';
            for iFile=1:numFiles
                plot(outData{iFile,1}, outData{iFile,2}(:,i) ...
                     ,'LineStyle',lStyle ...
                     ,'Marker',Markers{iFile} ...
                     ,'MarkerSize',3 ...
                     ,'DisplayName',[strrep(FASTfilesDesc{iFile},'\','\\') ' (' strrep(outData{iFile,3}{i},'_','\_') ')' ] ...
                     ,'Color',LineColors{iFile} ...
                     ,'LineWidth',LineWidthConst);
                hold on;      
                lStyle = ':';
            end
            ylabel({strrep(RefColumnTitles{i+1},'_','\_')     RefColumnUnits{i+1}},'FontSize',FntSz); %cell array = print on two lines            
        else
            f=figNo;
            figure(f);
            for iFile=1:numFiles
                plot(outData{iFile,1}, outData{iFile,2}(:,i) ...
                     ,'LineStyle','-' ...
                     ,'Marker',Markers{iFile} ...
                     ,'MarkerSize',3 ...
                     ,'DisplayName',[strrep(FASTfilesDesc{iFile},'\','\\') ' (' outData{iFile,3}{i} ', ' RefColumnUnits{i+1} ')'  ] ...
                     ,'LineWidth',LineWidthConst);
                hold on;      
            end   
            ylabel('FAST Channels','FontSize',FntSz);            
        end
        set(gca,'FontSize',FntSz-2,'gridlinestyle','-');
        xlabel([RefColumnTitles{1  } ' ' RefColumnUnits{1  }],'FontSize',FntSz); %string = print on one line
        title( strrep(titleText,'\','\\'),'FontSize',FntSz )    
        grid on;

        if numFiles > 1 && ShowLegend
            legend show %(FASTfilesDesc{:});
        end
% xlim([0,0.008])
% xlim([0,0.5])
        set(f,'Name',RefColumnTitles{i+1} ...
             ,'paperorientation','landscape' ...
             ,'paperposition',[0.25 0.25 10.5 8]);  
         
         
        if figNo < 0
             outFigName = [outFileRoot '_' num2str(i)];           
%            outFigName = ['AD_Results\AllPlots\' outFileRoot '_' num2str(i)];            
%             saveFigureAsPNG( f, outFigName  ) 
        end
        
    end


return
end

%%
function [] = plotPSDData( outData, FASTfilesDesc, Markers, LineColors, ...
                    RefColumnTitles, RefColumnUnits, titleText, ShowLegend, LineWidthConst, FntSz )
numCols  = size(outData{1,2},2) ;
numFiles = size(outData,1);

%% -----------------------------------------------------------
% Plot the psd from each file, with each channel in 
%      a separate figure:
% ------------------------------------------------------------
    for i = 1:numCols    
        f=figure;
        for iFile=1:numFiles
            
            n=length(outData{iFile,2}(:,i));
            if mod(n,2)==1
                n=n-1;
            end
            
            [ f1, Sf1 ] = getPSD( outData{iFile,2}(1:n,i), 1/(outData{iFile,1}(n) - outData{iFile,1}(1)) );
            loglog(f1, Sf1 ...
                 ,'LineStyle','-' ...
                 ,'Marker',Markers{iFile} ...
                 ,'MarkerSize',4 ...
                 ,'DisplayName',[strrep(FASTfilesDesc{iFile},'\','\\') ' (' strrep(outData{iFile,3}{i},'_','\_') ')' ] ...                 
                 ,'Color',LineColors{iFile} ...
                 ,'LineWidth',LineWidthConst);
            hold on;      
        end
        set(gca,'FontSize',FntSz-2,'gridlinestyle','-');
        ylabel({' PSD of ' strrep(RefColumnTitles{i+1},'_','\_')  RefColumnUnits{i+1}},'FontSize',FntSz); %cell array = print on two lines
        xlabel(['1/' RefColumnUnits{1  }],'FontSize',FntSz); %string = print on one line
        title( titleText,'FontSize',FntSz )    
        grid on;

        if numFiles > 1 && ShowLegend
            legend show %(FASTfilesDesc{:});
        end

        set(f,'Name',[RefColumnTitles{i+1} ' PSD']...
             ,'paperorientation','landscape' ...
             ,'paperposition',[0.25 0.25 10.5 8]);   
         
%          xlim([0,60]);
    end





return
end

%% possibly use this to make sure the channel names are the same....
%% ------------------------------------------------------------------------
function [Indx,err,ColToFind,scaleFact] = getColIndx( ColToFind, colNames, fileName )
    % ColToFind is originally the name of the reference channel
    % on exit, it is the channel name in the current file
    
    err = false;
    scaleFact = 1;
    Indx = find( strcmpi(ColToFind, colNames), 1, 'first' );
    
    if isempty(Indx) % let's try the negative of this column        
        
        [Indx,scaleFact,ColToFind] = getNewColIndx(ColToFind, colNames);
        
        if isempty(Indx)
            if strncmp(ColToFind,'-',1)
                ColToFind = ColToFind(2:end);
                scaleFact = -1;
            else
                ColToFind = strcat('-',ColToFind);
                scaleFact = -1;
            end        
            Indx = find( strcmpi(ColToFind, colNames), 1, 'first' );
            
            if isempty(Indx) % let's try the negative of this column              
                [Indx,scaleFact,ColToFind] = getNewColIndx(ColToFind, colNames);
            end            
        end
        
    end
    
    if isempty(Indx)
        disp(['Error: ' ColToFind ' not found in ' fileName ]);
        err = true;
    end
return
end


function [Indx,scaleFact,ColToFind] = getNewColIndx( ColToFind, colNames )
 
    scaleFact = 1;
%     [ColToFind,    scaleFact] = getAD14ChannelName(ColToFind);
%     [ColToFind,    scaleFact] = getFASTv7ChannelName(ColToFind);

    [newColNames, scaleFactM] = MooringColNames( ColToFind, colNames );
    [newColNames, scaleFactB] = BeamDynColNames( ColToFind, newColNames );    
    scaleFact = scaleFact*scaleFactM*scaleFactB;
    
%     disp( [colNames newColNames])

    Indx = find( strcmpi(ColToFind, newColNames), 1, 'first' );
    
    
    if isempty(Indx) % let's try the negative of this column        
        
        if strncmp(ColToFind,'-',1)
            ColToFind = ColToFind(2:end);
            scaleFact = -1*scaleFact;
        else
            ColToFind = strcat('-',ColToFind);
            scaleFact = -1*scaleFact;
        end        
        Indx = find( strcmpi(ColToFind, newColNames), 1, 'first' );
                    
    end
    
    if ~isempty(Indx)
        ColToFind = colNames{Indx};
    end
    
return
end

function [unitOut, scalefact] = ConvertUnits( unitIn )

    if strcmpi( unitIn, '(rad)' ) 
        unitOut = '(deg)';
        scalefact  = 180/pi;
    elseif strcmpi( unitIn, '(rad/s)' )
        unitOut = '(rpm)';
        scalefact  = 30/pi;
    elseif strcmpi( unitIn, '(rad/ss)' )
        unitOut = '(rpm/s)';
        scalefact  = 30/pi;
    elseif strcmpi( unitIn, '(N-m)' )
        unitOut = '(kN-m)';
        scalefact  = 1/1000;
    elseif strcmpi( unitIn, '(N)' )
        unitOut = '(kN)';
        scalefact  = 1/1000;
    else
        unitOut = unitIn;
        scalefact  = 1;        
    end    

end

function [colNames,scalefact] = BeamDynColNames( ColToFind, colNames )
    scalefact = 1;
    
    % from BD to ED names:
    if strcmpi(ColToFind, 'B1RootFxr') || ...
       strcmpi(ColToFind, 'B2RootFxr') || ...
       strcmpi(ColToFind, 'B3RootFxr')        
       colNames = strrep(colNames,'RootFxb1','B1RootFxr');
       colNames = strrep(colNames,'RootFxb2','B2RootFxr');
       colNames = strrep(colNames,'RootFxb3','B3RootFxr');
       scalefact = 1E3;      
    elseif strcmpi(ColToFind, 'B1RootFyr') || ...
           strcmpi(ColToFind, 'B2RootFyr') || ...
           strcmpi(ColToFind, 'B3RootFyr')        
       colNames = strrep(colNames,'RootFyb1','B1RootFyr');
       colNames = strrep(colNames,'RootFyb2','B2RootFyr');
       colNames = strrep(colNames,'RootFyb3','B3RootFyr');
       scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1RootFzr') || ...
           strcmpi(ColToFind, 'B2RootFzr') || ...
           strcmpi(ColToFind, 'B3RootFzr')        
       colNames = strrep(colNames,'RootFzb1','B1RootFzr');
       colNames = strrep(colNames,'RootFzb2','B2RootFzr');
       colNames = strrep(colNames,'RootFzb3','B3RootFzr');
       colNames = strrep(colNames,'RootFzc1','B1RootFzr');
       colNames = strrep(colNames,'RootFzc2','B2RootFzr');
       colNames = strrep(colNames,'RootFzc3','B3RootFzr');
       scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1RootMxr') || ...
           strcmpi(ColToFind, 'B2RootMxr') || ...
           strcmpi(ColToFind, 'B3RootMxr')        
       colNames = strrep(colNames,'RootMxb1','B1RootMxr');
       colNames = strrep(colNames,'RootMxb2','B2RootMxr');
       colNames = strrep(colNames,'RootMxb3','B3RootMxr');
       scalefact = 1E3;      
    elseif strcmpi(ColToFind, 'B1RootMyr') || ...
           strcmpi(ColToFind, 'B2RootMyr') || ...
           strcmpi(ColToFind, 'B3RootMyr')        
       colNames = strrep(colNames,'RootMyb1','B1RootMyr');
       colNames = strrep(colNames,'RootMyb2','B2RootMyr');
       colNames = strrep(colNames,'RootMyb3','B3RootMyr');
       scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1RootMzr') || ...
           strcmpi(ColToFind, 'B2RootMzr') || ...
           strcmpi(ColToFind, 'B3RootMzr')        
       colNames = strrep(colNames,'RootMzb1','B1RootMzr');
       colNames = strrep(colNames,'RootMzb2','B2RootMzr');
       colNames = strrep(colNames,'RootMzb3','B3RootMzr');
       colNames = strrep(colNames,'RootMzc1','B1RootMzr');
       colNames = strrep(colNames,'RootMzc2','B2RootMzr');
       colNames = strrep(colNames,'RootMzc3','B3RootMzr');
       scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1TipTDxr') || ...
           strcmpi(ColToFind, 'B2TipTDxr') || ...
           strcmpi(ColToFind, 'B3TipTDxr')        
       colNames = strrep(colNames,'TipDxb1','B1TipTDxr');
       colNames = strrep(colNames,'TipDxb2','B2TipTDxr');
       colNames = strrep(colNames,'TipDxb3','B3TipTDxr');
       scalefact = 1;
    elseif strcmpi(ColToFind, 'B1TipTDyr') || ...
           strcmpi(ColToFind, 'B2TipTDyr') || ...
           strcmpi(ColToFind, 'B3TipTDyr')        
       colNames = strrep(colNames,'TipDyb1','B1TipTDyr');
       colNames = strrep(colNames,'TipDyb2','B2TipTDyr');
       colNames = strrep(colNames,'TipDyb3','B3TipTDyr');
       scalefact = 1;
    elseif strcmpi(ColToFind, 'B1TipTDzr') || ...
           strcmpi(ColToFind, 'B2TipTDzr') || ...
           strcmpi(ColToFind, 'B3TipTDzr')        
       colNames = strrep(colNames,'TipDzb1','B1TipTDzr');
       colNames = strrep(colNames,'TipDzb2','B2TipTDzr');
       colNames = strrep(colNames,'TipDzb3','B3TipTDzr');
       colNames = strrep(colNames,'TipDzc1','B1TipTDzr');
       colNames = strrep(colNames,'TipDzc2','B2TipTDzr');
       colNames = strrep(colNames,'TipDzc3','B3TipTDzr');
       scalefact = 1;
    elseif strcmpi(ColToFind, 'B1TipTAXl') || ...
           strcmpi(ColToFind, 'B2TipTAXl') || ...
           strcmpi(ColToFind, 'B3TipTAXl')        
       colNames = strrep(colNames,'TipALxb1','B1TipTAXl');
       colNames = strrep(colNames,'TipALxb2','B2TipTAXl');
       colNames = strrep(colNames,'TipALxb3','B3TipTAXl');
       scalefact = 1; 
    elseif strcmpi(ColToFind, 'B1TipTAYl') || ...
           strcmpi(ColToFind, 'B2TipTAYl') || ...
           strcmpi(ColToFind, 'B3TipTAYl')        
       colNames = strrep(colNames,'TipALyb1','B1TipTAYl');
       colNames = strrep(colNames,'TipALyb2','B2TipTAYl');
       colNames = strrep(colNames,'TipALyb3','B3TipTAYl');
       scalefact = 1; 
    elseif strcmpi(ColToFind, 'B1TipTAZl') || ...
           strcmpi(ColToFind, 'B2TipTAZl') || ...
           strcmpi(ColToFind, 'B3TipTAZl')        
       colNames = strrep(colNames,'TipALzb1','B1TipTAZl');
       colNames = strrep(colNames,'TipALzb2','B2TipTAZl');
       colNames = strrep(colNames,'TipALzb3','B3TipTAZl');
       scalefact = 1; 
    elseif strcmpi(ColToFind, 'B1N2Fxl') || ...
           strcmpi(ColToFind, 'B2N2Fxl') || ...
           strcmpi(ColToFind, 'B3N2Fxl') 
        colNames = strrep(colNames,'Spn1FLxb1','B1N2Fxl');
        colNames = strrep(colNames,'Spn1FLxb2','B2N2Fxl');
        colNames = strrep(colNames,'Spn1FLxb3','B3N2Fxl');
        scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1N2Fyl') || ...
           strcmpi(ColToFind, 'B2N2Fyl') || ...
           strcmpi(ColToFind, 'B3N2Fyl')  
        colNames = strrep(colNames,'Spn1FLyb1','B1N2Fyl');  
        colNames = strrep(colNames,'Spn1FLyb2','B2N2Fyl');  
        colNames = strrep(colNames,'Spn1FLyb3','B3N2Fyl');  
        scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1N2Fzl') || ...
           strcmpi(ColToFind, 'B2N2Fzl') || ...
           strcmpi(ColToFind, 'B3N2Fzl') 
        colNames = strrep(colNames,'Spn1FLzb1','B1N2Fzl');
        colNames = strrep(colNames,'Spn1FLzb2','B2N2Fzl');
        colNames = strrep(colNames,'Spn1FLzb3','B3N2Fzl');
        scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1N2Mxl') || ...
           strcmpi(ColToFind, 'B2N2Mxl') || ...
           strcmpi(ColToFind, 'B3N2Mxl') 
        colNames = strrep(colNames,'Spn1MLxb1','B1N2Mxl');
        colNames = strrep(colNames,'Spn1MLxb2','B2N2Mxl');
        colNames = strrep(colNames,'Spn1MLxb3','B3N2Mxl');
        scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1N2Myl') || ...
           strcmpi(ColToFind, 'B2N2Myl') || ...
           strcmpi(ColToFind, 'B3N2Myl')  
        colNames = strrep(colNames,'Spn1MLyb1','B1N2Myl');
        colNames = strrep(colNames,'Spn1MLyb2','B2N2Myl');
        colNames = strrep(colNames,'Spn1MLyb3','B3N2Myl');
        scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1N2Mzl') || ...
           strcmpi(ColToFind, 'B2N2Mzl') || ...
           strcmpi(ColToFind, 'B3N2Mzl') 
        colNames = strrep(colNames,'Spn1MLzb1','B1N2Mzl');
        colNames = strrep(colNames,'Spn1MLzb2','B2N2Mzl');
        colNames = strrep(colNames,'Spn1MLzb3','B3N2Mzl');
        scalefact = 1E3;

    elseif strcmpi(ColToFind, 'B1N3Fxl') || ...
           strcmpi(ColToFind, 'B2N3Fxl') || ...
           strcmpi(ColToFind, 'B3N3Fxl') 
        colNames = strrep(colNames,'Spn2FLxb1','B1N3Fxl');
        colNames = strrep(colNames,'Spn2FLxb2','B2N3Fxl');
        colNames = strrep(colNames,'Spn2FLxb3','B3N3Fxl');
        scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1N3Fyl') || ...
           strcmpi(ColToFind, 'B2N3Fyl') || ...
           strcmpi(ColToFind, 'B3N3Fyl') 
        colNames = strrep(colNames,'Spn2FLyb1','B1N3Fyl');
        colNames = strrep(colNames,'Spn2FLyb2','B2N3Fyl');
        colNames = strrep(colNames,'Spn2FLyb3','B3N3Fyl');
        scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1N3Fzl') || ...
           strcmpi(ColToFind, 'B2N3Fzl') || ...
           strcmpi(ColToFind, 'B3N3Fzl') 
        colNames = strrep(colNames,'Spn2FLzb1','B1N3Fzl');
        colNames = strrep(colNames,'Spn2FLzb2','B2N3Fzl');
        colNames = strrep(colNames,'Spn2FLzb3','B3N3Fzl');
        scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1N3Mxl') || ...
           strcmpi(ColToFind, 'B2N3Mxl') || ...
           strcmpi(ColToFind, 'B3N3Mxl') 
        colNames = strrep(colNames,'Spn2MLxb1','B1N3Mxl');
        colNames = strrep(colNames,'Spn2MLxb2','B2N3Mxl');
        colNames = strrep(colNames,'Spn2MLxb3','B3N3Mxl');
        scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1N3Myl') || ...
           strcmpi(ColToFind, 'B2N3Myl') || ...
           strcmpi(ColToFind, 'B3N3Myl') 
        colNames = strrep(colNames,'Spn2MLyb1','B1N3Myl');
        colNames = strrep(colNames,'Spn2MLyb2','B2N3Myl');
        colNames = strrep(colNames,'Spn2MLyb3','B3N3Myl');
        scalefact = 1E3;
    elseif strcmpi(ColToFind, 'B1N3Mzl') || ...
           strcmpi(ColToFind, 'B2N3Mzl') || ...
           strcmpi(ColToFind, 'B3N3Mzl') 
        colNames = strrep(colNames,'Spn2MLzb1','B1N3Mzl');
        colNames = strrep(colNames,'Spn2MLzb2','B2N3Mzl');
        colNames = strrep(colNames,'Spn2MLzb3','B3N3Mzl');
        scalefact = 1E3;
        
    end
            
end


  

function [colNames,scalefact] = MooringColNames( ColToFind, colNames )
    scalefact = 1;
    
    if strncmpi(ColToFind, 'TFair',5) || strncmpi(ColToFind, 'TAnch',5) || ...
       strncmpi(ColToFind,'-TFair',6) || strncmpi(ColToFind,'-TAnch',6) || ...
       strncmpi(ColToFind,'MTFair',6) || strncmpi(ColToFind,'MTAnch',6)
          
%        colNames = strrep(colNames,'TFair','T');
%        colNames = strrep(colNames,'TAnch','T_a');
       colNames = strrep(colNames,'FAIRTEN','TFair[');
       colNames = strrep(colNames,'ANCHTEN','TAnch['); %Bjj: watch out for case
       for i=1:length(colNames)
           colNames{i} = strcat(colNames{i},']');
       end       
       scalefact = 1E-3;
       
    elseif strncmpi(ColToFind, 'FairTen',7) || strncmpi(ColToFind, 'AnchTen',7) || ...
           strncmpi(ColToFind,'-FairTen',8) || strncmpi(ColToFind,'-AnchTen',8) || ...
           strncmpi(ColToFind,'MFairTen',8) || strncmpi(ColToFind,'MAnchTen',8)
           
%        colNames = strrep(colNames,'TFair[','FairTen');
%        colNames = strrep(colNames,'TAnch[','AnchTen');   
%       scalefact = 1E3;
        colNames = strrep(colNames,'T[','FairTen');
        colNames = strrep(colNames,'T_a[','AnchTen');   
        scalefact = 1;
       colNames = strrep(colNames,']',''); 
       
       
    elseif strncmpi(ColToFind, 'T[',2) || strncmpi(ColToFind, 'T_a[',4) || ...
           strncmpi(ColToFind,'-T[',3) || strncmpi(ColToFind,'-T_a[',5) || ...
           strncmpi(ColToFind,'MT[',3) || strncmpi(ColToFind,'MT_a[',5)
           
       colNames = strrep(colNames,'FAIRTEN','T[');
       colNames = strrep(colNames,'ANCHTEN','T_a[');   
       for i=1:length(colNames)
           colNames{i} = strcat(colNames{i},']');
       end       
       scalefact = 1;
       
    end
    
end


        
function [ChannelName_old,scaleFact] = getFASTv7ChannelName(ChannelName)

    scaleFact = 1.0;

    if strcmpi(ChannelName,'Wave1Elev')
        ChannelName_old = 'WaveElev';                        
    elseif strcmpi(ChannelName,'TwHt1TPxi')
        ChannelName_old = 'TTDspFA';            
    elseif strcmpi(ChannelName,'uWind')
        ChannelName_old = 'Wind1VelX';            
    elseif strcmpi(ChannelName,'vWind')
        ChannelName_old = 'Wind1VelY';            
    elseif strcmpi(ChannelName,'wWind')
        ChannelName_old = 'Wind1VelZ';            
    elseif strcmpi(ChannelName,'RotCp')
        ChannelName_old = 'RtAeroCp';            
        scaleFact = 1;
    elseif strcmpi(ChannelName,'RotCq')
        ChannelName_old = 'RtAeroCq';            
        scaleFact = 1;
    elseif strcmpi(ChannelName,'TSR')
        ChannelName_old = 'RtTSR';            
    elseif strcmpi(ChannelName,'RotCt')
        ChannelName_old = 'RtAeroCt';            
        scaleFact = 1;
    elseif strcmpi(ChannelName,'TwHt1TPyi')
        ChannelName_old = 'TTDspSS';                    
    elseif strcmpi(ChannelName,'ReactFXss')
        ChannelName_old = '-TwrBsFxt';
        scaleFact = 1000;
    elseif strcmpi(ChannelName,'ReactFYss')
        ChannelName_old = '-TwrBsFyt';
        scaleFact = 1000;
    elseif strcmpi(ChannelName,'ReactFZss')
        ChannelName_old = '-TwrBsFzt';
        scaleFact = 1000;
    elseif strcmpi(ChannelName,'ReactMXss')
        ChannelName_old = '-TwrBsMxt';
        scaleFact = 1000;
    elseif strcmpi(ChannelName,'ReactMYss')
        ChannelName_old = '-TwrBsMyt';
        scaleFact = 1000;
    elseif strcmpi(ChannelName,'ReactMZss')
        ChannelName_old = '-TwrBsMzt';       
        scaleFact = 1000;                                                                                  
    else
        ChannelName_old = strrep(ChannelName,'M1N1V','Wave1V');
        ChannelName_old = strrep(ChannelName_old,'M1N1A','Wave1A');
    end
                                
    return 

end

function [ChannelName_old,scaleFact] = getOldFASTv8ChannelName(ChannelName)

    scaleFact = 1.0;

    if strcmpi(ChannelName,'IntfFXss')
        ChannelName_old = 'IntfXss';            
    elseif strcmpi(ChannelName,'IntfFYss')
        ChannelName_old = 'IntfYss';            
    elseif strcmpi(ChannelName,'IntfFZss')
        ChannelName_old = 'IntfZss';            
    elseif strcmpi(ChannelName,'ReactFXss')
        ChannelName_old = 'ReactXss';            
    elseif strcmpi(ChannelName,'ReactFYss')
        ChannelName_old = 'ReactYss';            
    elseif strcmpi(ChannelName,'ReactFZss')
        ChannelName_old = 'ReactZss'; 
    else
        ChannelName_old = ChannelName;
    end
                                
    return 

end


function [ChannelName_new,scaleFact] = getAD14ChannelName(ChannelName)
       
    scaleFact = 1.0;
    ChannelName_new = ChannelName;
        
    % there isn't an easy way to get this node (i'd need to read the
    % AD15 input file to do so)
        
%     AD14node = '06';  %for test01
%    dr = 1.2573;
%     AD14node = '10'; %for test 10
%     dr = 2.3490000E-01;
    
%     AD14node = '07'; %for test 20
%     dr = 4.1000000E+00;
%     
%     AD14node = '07'; %for test 12
%     dr=2.2166700E+00;
    
%     TestNodes = [3,9,11,13,14,15,16,17,18];
%     TestNodes = [3,9,11,13,14,15,16,17,18];
%     drNodes = [2.7333000E+00
% 2.7333000E+00
% 2.7333000E+00
% 4.1000000E+00
% 4.1000000E+00
% 4.1000000E+00
% 4.1000000E+00
% 4.1000000E+00
% 4.1000000E+00
% 4.1000000E+00
% 4.1000000E+00
% 4.1000000E+00
% 4.1000000E+00
% 4.1000000E+00
% 2.7333000E+00
% 2.7333000E+00
% 2.7333000E+00];
%SWRT
    TestNodes = [2, 8, 16];
    drNodes = [1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01
1.7310000E-01];




    
    if length(ChannelName) < 5 
        return;
    end

    if strcmpi(ChannelName(1:3),'B1N')
        
        n=eval(ChannelName(4));
        AD14n = TestNodes(n)-1;
        dr = drNodes(AD14n);
        AD14node = sprintf('%02.0f',AD14n);
                    
        B1N2_channel = ChannelName(5:end);

        if strcmpi(B1N2_channel,'Alpha')
            ChannelName_new = [ 'Alpha' AD14node];
        elseif strcmpi(B1N2_channel,'DynP')
            ChannelName_new = [ 'DynPres' AD14node];
        elseif strcmpi(B1N2_channel,'Cl')
            ChannelName_new = [ 'Clift' AD14node];
        elseif strcmpi(B1N2_channel,'Cd')
            ChannelName_new = [ 'Cdrag' AD14node];
        elseif strcmpi(B1N2_channel,'Cn')
            ChannelName_new = [ 'CNorm' AD14node];
        elseif strcmpi(B1N2_channel,'Ct')
            ChannelName_new = [ 'CTang' AD14node];
        elseif strcmpi(B1N2_channel,'Cm')
            ChannelName_new = [ 'CMomt' AD14node];
        elseif strcmpi(B1N2_channel,'Theta')
            ChannelName_new = [ 'Pitch' AD14node];
        elseif strcmpi(B1N2_channel,'AxInd')
            ChannelName_new = [ 'AxInd' AD14node];
        elseif strcmpi(B1N2_channel,'TnInd')
            ChannelName_new = [ 'TanInd' AD14node];
         elseif strcmpi(B1N2_channel,'Fx')
            ChannelName_new = [ 'ForcN' AD14node];
            scaleFact = 1/dr;
         elseif strcmpi(B1N2_channel,'Fy')
            ChannelName_new = [ 'ForcT' AD14node];
            scaleFact = 1/dr;
         elseif strcmpi(B1N2_channel,'Mm')
            ChannelName_new = [ 'Pmomt' AD14node];
            scaleFact = 1/dr;
         elseif strcmpi(B1N2_channel,'Re')
            ChannelName_new = [ 'ReNum' AD14node];
         elseif strcmpi(B1N2_channel,'Vx')
            ChannelName_new = [ 'Vx' AD14node];
         elseif strcmpi(B1N2_channel,'Vy')
            ChannelName_new = [ 'Vy' AD14node];
        end

    end

    return

end


function [ChannelName_new,scaleFact] = getNewChannelName(ChannelName)
        
        scaleFact = 1.0;
        ChannelName_new = ChannelName;
        
        
%         if strcmpi(ChannelName,'WaveElev')
%             ChannelName_new = 'Wave1Elev';
%         elseif strncmpi( ChannelName,'Fair',4 ) && strcmpi( ChannelName((end-2):end), 'Ten') %TFair[i] = FairiTen
%             ChannelName_new = strrep( strrep( ChannelName,'Fair','TFair['),'Ten',']');
%         elseif strncmpi( ChannelName,'Anch',4 ) && strcmpi( ChannelName((end-2):end), 'Ten') %TAnch[i] = AnchiTen
%             ChannelName_new = strrep( strrep( ChannelName,'Anch','TAnch['),'Ten',']');
%         elseif strncmpi( ChannelName,'Anch',4 ) && strcmpi( ChannelName((end-2):end), 'Ten') %TAnch[i] = AnchiTen
%             ChannelName_new = strrep( strrep( ChannelName,'Anch','TAnch['),'Ten',']');
%         elseif strcmpi(ChannelName,'IntfXss')
%             ChannelName_new = 'IntfFXss';            
%         elseif strcmpi(ChannelName,'IntfYss')
%             ChannelName_new = 'IntfFYss';            
%         elseif strcmpi(ChannelName,'IntfZss')
%             ChannelName_new = 'IntfFZss';            
%         elseif strcmpi(ChannelName,'ReactXss')
%             ChannelName_new = 'ReactFXss';            
%         elseif strcmpi(ChannelName,'ReactYss')
%             ChannelName_new = 'ReactFYss';            
%         elseif strcmpi(ChannelName,'ReactZss')
%             ChannelName_new = 'ReactFZss';            
%             
%             
%         elseif strcmpi(ChannelName,'TTDspFA')
%             ChannelName_new = 'TwHt1TPxi';            
%         elseif strcmpi(ChannelName,'TTDspSS')
%             ChannelName_new = 'TwHt1TPyi';                    
%         elseif strcmpi(ChannelName,'-TwrBsFxt')
%             ChannelName_new = 'ReactFXss';
%             scaleFact = 1000;
%         elseif strcmpi(ChannelName,'-TwrBsFyt')
%             ChannelName_new = 'ReactFYss';
%             scaleFact = 1000;
%         elseif strcmpi(ChannelName,'-TwrBsFzt')
%             ChannelName_new = 'ReactFZss';
%             scaleFact = 1000;
%         elseif strcmpi(ChannelName,'-TwrBsMxt')
%             ChannelName_new = 'ReactMXss';
%             scaleFact = 1000;
%         elseif strcmpi(ChannelName,'-TwrBsMyt')
%             ChannelName_new = 'ReactMYss';
%             scaleFact = 1000;
%         elseif strcmpi(ChannelName,'-TwrBsMzt')
%             ChannelName_new = 'ReactMZss';       
%             scaleFact = 1000;
%             
%             
%                                                           
%         else
%             ChannelName_new = strrep(ChannelName,'Wave1V','M1N1V');
%             ChannelName_new = strrep(ChannelName_new,'Wave1A','M1N1A');
%             %ChannelName_new = ChannelName_new;              
%         end
                        
        
return     
end 

%% ------------------------------------------------------------------------
