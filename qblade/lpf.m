Fs=20;
d=designfilt('lowpassfir','FilterOrder',10,'cutoffFrequency',1e-10,'SampleRate',Fs)

figure
freqz(d.Coefficients, 1, 2^14, Fs)

d = designfilt("lowpassfir", ...
    PassbandFrequency=1e-3,StopbandFrequency=1e-3, ...
    PassbandRipple=1,StopbandAttenuation=60, ...
    DesignMethod="equiripple");
y = filtfilt(d,x);

load('')

d1 = designfilt("lowpassiir",FilterOrder=45, ...
    HalfPowerFrequency=0.15,DesignMethod="butter");
y = filtfilt(d1,rotSpeedFull);

figure
plot(rotSpeedFull)
hold on
plot(y)