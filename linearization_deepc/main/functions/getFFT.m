function [f,uAmpSpectrum,uPSD] = getFFT(Ts,u)
%getFFT(Ts,u) A funtion used to get the spectrum of the input signal u
%sampled with sampling period Ts.
%
% Input arguments
%----------------
% Ts    : Sampling time (in seconds).
% u     : Input signal.
%
% Output arguments: 
%------------------
% f                : frequencies vector.
% uAmpSpectrum     : amplitude of spectrum.
% uPSD             : PSD estimate.
%==========================================================================

N = max(size(u));
fs = 1/Ts; % Sampling frequency (in Hz)
f = fs/2*linspace(0,1,N/2+1); % frequencies vector
uFFT = fft(u);
uAmpSpectrum = abs(uFFT/N); % Normalize and extract amplitude only
uAmpSpectrum = uAmpSpectrum(1:N/2+1); % Truncate negative frequencies
uAmpSpectrum(2:end-1) = 2*uAmpSpectrum(2:end-1); % Compensate for truncated frequencies, DC component and Nyquist frequency do not occur twice

% figure
% plot(f,uAmpSpectrum)
% title('Amplitude Spectrum')
% xlabel('Frequency (in Hz)')
% ylabel('Amplitude')
% xlim([0 1])
% grid on
% set(gcf,'Color','White')

% PSD estimate
uPSD = pow2db(2.*(abs(uFFT(1:N/2+1)).^2/(N*fs)));

% figure
% periodogram(u,rectwin(N),N,fs)
% % plot(f,pow2db(uPSD))
% title("Power Spectrum Density Estimate")
% xlabel("Frequency (Hz)")
% ylabel("Power/Frequency (dB/Hz)")
% grid on
% set(gcf,'Color','White')

end