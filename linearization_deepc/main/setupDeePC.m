function [SysData, controlParams, figHandle] = setupDeePC(G, Ts, RandomSignalInfo, ctrlInputIdx, previewFlag, prevInputIdx, measNoiseSTD, N, f, p, plotIO)
% setupDeePC(G,Ts, contrInpName, RandomSignalInfo, scalingFactors)
% Function used to obtain signals for open loop training of DeePC and
% training of DeePC.
%
% Input arguments:
%-----------------
% G             : System
% Ts            : Sampling time
% simLen        : Length of produced signal
% RandomSignalInfo : Struct with fields :   Type: type of random signal 'constZero' | 'rbs' | 'rgs' | 'prbs' | 'sine'
%                                           Band: Frequency range of generated signal
%                                           Range: Range of random signal
%                                           ScalingFactor: Factor that the random signal is multiplied with after generation
% ctrlInputIdx  : Index in RandomSignalInfo of the controllable signal
% prevInputIdx  : Indices in RandomSignalInfo of the signals with preview
% N             : Lenght of data set stored for closed loop control
% f             : Future window length
% p             : Past window length
% plotIO        : Boolean whether to plot or not
%
% Output arguments:
%------------------
% SysData       : Data obtained from the system
% controlParams : Data for control schemes
% figHandles    : Figure handles if plotIO == true, empty array otherwise
%==========================================================================

% Discretize system
G_d = c2d(G,Ts);
SysData.discreteSS = G_d;

% Time vector
t = (0:N*2)*Ts;
nt = length(t);

% Generate signal per system input
nu = length(RandomSignalInfo);
u = zeros(nt, nu);
for idx = 1:length(RandomSignalInfo)
    if ~strcmp(RandomSignalInfo(idx).Type, 'constZero')
        u(:,idx) = idinput(nt,...
            RandomSignalInfo(idx).Type,...
            RandomSignalInfo(idx).Band,...
            RandomSignalInfo(idx).Range)...
            .*RandomSignalInfo(idx).ScalingFactor;
    end
    SysData.OL.input(idx).inputName = RandomSignalInfo(idx).inputName;
    SysData.OL.input(idx).signal = u(:,idx);
end

% Simulate system and get open-loop I/O data set
y = lsim(G_d, u, t); % zero initial condition, around linearization point

% Add measurement noise on output
if measNoiseSTD ~= 0
    y = y + measNoiseSTD.*randn(size(y));
end
SysData.OL.output = y;

if plotIO
    figHandle = figure;
    subplot(1,2,1)
    plot(t,y)
    grid on
    xlabel('Time (in s)')
    ylabel('Rotor speed, scaled (-)')
    title('OL response to PE input - \omega around linearization point')
    set(gcf,'Color','White')
    
    subplot(1,2,2)
    yyaxis left
    plot(t,u(:,prevInputIdx))
    ylabel('Disturbance')
    grid on
    yyaxis right
    plot(t,u(:,ctrlInputIdx))
    ylabel('Input')
    grid on
    xlabel('Time (in s)')
    title('OL input & disturbances')
    set(gcf,'Color','White')
    
end

%% DeePC parameters
Nbar = N-p-f+1;
i = 1;

controlParams.N = N;
controlParams.p = p;
controlParams.f = f;

%% Construct data matrices
inputData = zeros(nt, length(ctrlInputIdx));
for i=1:length(ctrlInputIdx); inputData(:,i) = SysData.OL.input(ctrlInputIdx(i)).signal; end
% Past data
SysData.Up = constructHankelMat(inputData,i,p,Nbar);
SysData.Yp = constructHankelMat(SysData.OL.output,i,p,Nbar);

% Future data
SysData.Uf = constructHankelMat(inputData,i+p,f,Nbar);
SysData.Yf = constructHankelMat(SysData.OL.output,i+p,f,Nbar);

disturbMat = u(:, prevInputIdx);

if previewFlag == 1
    SysData.Wp = constructHankelMat(disturbMat,i,p,Nbar); % past data
    SysData.Wf = constructHankelMat(disturbMat,i+p,f,Nbar); % future data
else
    SysData.Wp = [];
    SysData.Wf = [];
end


end



