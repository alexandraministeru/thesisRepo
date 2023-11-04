load('inputData\waveForces.mat');
load('inputData\waveLTIsys_reduced_rotSpeed.mat')

scaling = 1;

if scaling == 1
    % Scaling factors
    uhat_max = 10*(pi/180); % Maximum expected input (rad)
    v0hat_max = 0.1*16; % Maximum expected wind disturbance (m/s)
    FsurgHat_max = max(F_surge); % Maximum expected wave surge force (N)
    MpitchHat_max = max(M_pitch); % Maximum expected wave pitch moment (Nm)
    ehat_max = 0.1* 1173.69; % Maximum expected generator speed error (10% around linearization OP) (rpm)

    % Du = diag([uhat_max v0hat_max FsurgHat_max MpitchHat_max]);
    Du = diag([uhat_max v0hat_max MpitchHat_max]);
    Dy = ehat_max;

    % Scaled transfer functions
    G = Dy\LTIsys_reduced*Du;
end

fWaveMin = 0.05; % Hz
fWaveMax = 0.5; % Hz

figure
bPlot1 = bodeplot(LTIsys_reduced(1,1));
setoptions(bPlot1,'FreqUnits','Hz')
grid on
set(gcf,'Color','White')

figure
bPlot2 = bodeplot(-LTIsys_reduced(1,2));
setoptions(bPlot2,'FreqUnits','Hz')
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
grid on
set(gcf,'Color','White')

figure
bPlot3 = bodeplot(-LTIsys_reduced(1,3));
setoptions(bPlot3,'FreqUnits','Hz')
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
grid on
set(gcf,'Color','White')

%% del
figure
bPlot1 = bodeplot(LTIsys_reduced(2,1));
setoptions(bPlot1,'FreqUnits','Hz')
grid on
set(gcf,'Color','White')

figure
bPlot2 = bodeplot(-LTIsys_reduced(2,2));
setoptions(bPlot2,'FreqUnits','Hz')
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
grid on
set(gcf,'Color','White')

figure
bPlot3 = bodeplot(-LTIsys_reduced(2,3));
setoptions(bPlot3,'FreqUnits','Hz')
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
grid on
set(gcf,'Color','White')

figure
step(LTIsys_reduced(2,1))
[v,t] = step(LTIsys_reduced(2,1));
sv = stepinfo(v,t)

figure
step(-LTIsys_reduced(2,2))
[fs,t] = step(-LTIsys_reduced(2,2));
sfs = stepinfo(fs,t)

figure
step(LTIsys_reduced(2,3))
[mp,t] = step(LTIsys_reduced(2,3));
smp = stepinfo(mp,t)

%% Step response
figure
step(-LTIsys_reduced(1,1))
[v,t] = step(-LTIsys_reduced(1,1));
sv = stepinfo(v,t)

figure
step(LTIsys_reduced(1,2))
[fs,t] = step(LTIsys_reduced(1,2));
sfs = stepinfo(fs,t)

figure
step(-LTIsys_reduced(1,3))
[mp,t] = step(-LTIsys_reduced(1,3));
smp = stepinfo(mp,t)
%%
Ts = 0.05;

% Turbulent wind
load('inputData\turbWind_16mps_long.mat')
v = turbWind;
[f_v,fft_v,psd_v] = getFFT(Ts,v);
tsim = 0:Ts:Ts*(numel(v)-1);

figure
sp1 = subplot(2,1,1)
plot(f_v(find(f_v==0):end),fft_v(find(f_v==0):end))
xlabel('Frequency (in Hz)')
ylabel('Amplitude')
title('FFT')
xlim([0 0.2])
grid on
sp2 = subplot(2,1,2)
plot(f_v(find(f_v==0):end),psd_v(find(f_v==0):end))
xlabel('Frequency (in Hz)')
ylabel('Power/Frequency (dB/Hz)')
title('PSD')
xlim([0 0.2])
grid on
linkaxes([sp1,sp2],'x');
sgtitle('Wind')
set(gcf,'Color','White')


[f_Fsg,fft_Fsg,psd_Fsg] = getFFT(Ts,F_surge);
[f_Mp,fft_Mp,psd_Mp] = getFFT(Ts,M_pitch);
tsim = 0:Ts:Ts*(numel(F_surge)-1);

fSurgePtfm = 0.051/(2*pi); %Hz
fPitchPtfm = 0.215/(2*pi); %Hz

figure
sp1 = subplot(2,1,1)
plot(f_Fsg(find(f_Fsg==0):end),fft_Fsg(find(f_Fsg==0):end))
xlabel('Frequency (in Hz)')
ylabel('Amplitude')
title('FFT')
% xline(fSurgePtfm,'k--','Platform surge')
% xline(fPitchPtfm,'k--','Platform pitch')
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
xlim([0 1])
grid on
sp2 = subplot(2,1,2)
plot(f_Fsg(find(f_Fsg==0):end),psd_Fsg(find(f_Fsg==0):end))
xlabel('Frequency (in Hz)')
ylabel('Power/Frequency (dB/Hz)')
title('PSD')
% xline(fSurgePtfm,'k--','Platform surge')
% xline(fPitchPtfm,'k--','Platform pitch')
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
xlim([0 1])
grid on
linkaxes([sp1,sp2],'x');
sgtitle('Surge force')
set(gcf,'Color','White')


figure
sp1 = subplot(2,1,1)
plot(f_Mp(find(f_Mp==0):end),fft_Mp(find(f_Mp==0):end))
xlabel('Frequency (in Hz)')
ylabel('Amplitude')
title('FFT')
% xline(fSurgePtfm,'k--','Platform surge')
% xline(fPitchPtfm,'k--','Platform pitch')
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
xlim([0 1])
grid on
sp2 = subplot(2,1,2)
plot(f_Mp(find(f_Mp==0):end),psd_Mp(find(f_Mp==0):end))
xlabel('Frequency (in Hz)')
ylabel('Power/Frequency (dB/Hz)')
title('PSD')
% xline(fSurgePtfm,'k--','Platform surge')
% xline(fPitchPtfm,'k--','Platform pitch')
% xline(fWaveMin,'k:','WaveMin')
% xline(fWaveMax,'k:','WaveMax')
xlim([0 1])
grid on
linkaxes([sp1,sp2],'x');
sgtitle('Pitch moment')
set(gcf,'Color','White')


