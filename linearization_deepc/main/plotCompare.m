file1 = 'outputData\Fsg_V_NP';
file2 = 'outputData\Fsg_V_WP';

% % Update data structure for old data files
% inName = '\theta_c (in deg)';
% outName = 'Generator speed (in rpm)';
% scaledFlag = 0;
% controlParams.lbu = -10;
% controlParams.ubu = 10;
% controlParams.duDeg = 8;
% controlParams.Q = Q;
% description = '';
% save(file1,'inName','outName','scaledFlag','controlParams','description', ...
% "-append")

plotComp(file1,file2)

% Required: 
% inName: name and measurement unit of control input
% outName: name and measurement unit of controlled output
% tsim: simulation time vector
% ref: reference vector
% Ts: sampling period
% controlParams: p,f,Q,R
% out: controlled output sequence
% uSeq: control input sequence
% descriprion: short description mentioning preview content if preview is
%              enabled
% scaledFlag: flag indicating whether the system was scaled or not
% uhat_max, ehat_max: only necessary if scaledFlag == 1, indicate scaling
%                     values for signal reconstruction to real life values


function [] = plotComp(varargin)
% plotComp(varargin) Function used to compare control output, control 
% input and control input rate of an arbitrary number of simulations.
%
% Input arguments:
%-----------------
% varargin : file names containing data to be compared. If the type of
% comparison is No preview VS Preview or Not scaled VS Scaled, the files
% must be provided in this order.
%
% Output arguments:
%------------------
% -
%==========================================================================

nFiles = length(varargin);
if nFiles == 0
    error('No files provided.')
end

data = cell(nFiles,1);
for idxFile = 1:nFiles
    data{idxFile} = load(varargin{idxFile});
end

kFinal = data{1}.kFinal;
Ts = data{1}.Ts;
f = data{1}.f;
legendStr = cell(nFiles,1);

contentAns = input(['Plot comparison type:\n 1 - No preview VS Preview\n ...' ...
    '2 - Different preview components\n ...' ...
    '3 - Tuning\n ...' ...
    '4 - Not scaled VS Scaled: ']);
if contentAns == 3
    tunedVar = input('Which variable? Q/ R/ p/ f: ','s');
end

% Build legend
for idxOut = 1:nFiles
    switch contentAns
        case 2
            legendStr{idxOut} = data{idxOut}.description;
        case 3 
            switch tunedVar
                case {'Q','R'}
                    legendStr{idxOut} = [tunedVar '=' num2str(eval(['data{' num2str(idxOut) '}.controlParams.' tunedVar '(1,1)']))];
                otherwise
                    legendStr{idxOut} = [tunedVar '=' num2str(eval(['data{' num2str(idxOut) '}.' tunedVar]))];                    
            end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for idxOut = 1:nFiles
    % Reconstruct data if it was scaled
    if data{idxOut}.scaledFlag == 1
        data{idxOut}.out = data{idxOut}.out.*data{idxOut}.ehat_max;
    end
    plot(data{idxOut}.tsim,data{idxOut}.out);    
    hold on
end
xlabel('Time (in s)')
ylabel(data{1}.outName,Interpreter='latex')
title('DeePC for reference tracking')
grid on
hold on
plot(data{1}.tsim,data{1}.ref(1:kFinal)) % reference
xline(Ts*f,'k--','Future window size')
% stepIdxs = find(data{1}.ref);
% xline(Ts*stepIdxs(1),'k--','Reference step')
switch contentAns
    case 1
        legend('No preview','With preview','Reference','Location','SouthEast')
    case {2,3}
        legend(legendStr,'Location','SouthEast')
    case 4
        legend('Not scaled','Scaled','Reference','Location','SouthEast')
end
set(gcf,'Color','White')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for idxOut = 1:nFiles
    % Reconstruct data if it was scaled
    if data{idxOut}.scaledFlag == 1
        data{idxOut}.uSeq = rad2deg(data{idxOut}.uSeq.*data{idxOut}.uhat_max);
    end
    plot(data{idxOut}.tsim,rad2deg(data{idxOut}.uSeq));    
    hold on
end
xlabel('Time (in s)')
ylabel(data{1}.inName,Interpreter='latex')
ylim([-15 15])
yline(data{1}.controlParams.lbu,'r--','LineWidth',1)
yline(data{1}.controlParams.ubu,'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
% stepIdxs = find(data{1}.ref);
% xline(Ts*stepIdxs(1),'k--','Reference step')
title('Control input')
switch contentAns
    case 1
        legend('No preview','With preview','Reference','Location','SouthEast')
    case {2,3}
        legend(legendStr,'Location','SouthEast')
    case 4
        legend('Not scaled','Scaled','Reference','Location','SouthEast')
end
grid on
set(gcf,'Color','White')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% Control input rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for idxOut = 1:nFiles
    % Reconstruct data if it was scaled
    if data{idxOut}.scaledFlag == 1
        plot(data{idxOut}.tsim(1:end-1),diff(rad2deg(data{idxOut}.uSeq))*(1/Ts).*data{idxOut}.uhat_max); 
    else
        plot(data{idxOut}.tsim(1:end-1),diff(rad2deg(data{idxOut}.uSeq))*(1/Ts)); 
    end
    hold on
end
xlabel('Time (in s)')
ylabel('\theta_c rate (in deg/s)',Interpreter='latex')
yline(data{1}.controlParams.duDeg,'r--','LineWidth',1)
yline(-data{1}.controlParams.duDeg,'r--','LineWidth',1)
xline(Ts*f,'k--','Future window size')
ylim([-15 15])
title('Control input rate')
switch contentAns
    case 1
        legend('No preview','With preview','Reference','Location','SouthEast')
    case {2,3}
        legend(legendStr,'Location','SouthEast')
    case 4
        legend('Not scaled','Scaled','Reference','Location','SouthEast')
end
grid on
set(gcf,'Color','White')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

