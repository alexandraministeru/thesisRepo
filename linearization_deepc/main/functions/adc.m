function ADC = adc(uSeq,Ts,maxPitchRate)
pitchingRate = diff(uSeq).*(1/Ts);
timeWindow = 0:Ts:(length(uSeq)-1)*Ts;
n = numel(timeWindow);
ADC = 0;

for idxADC=1:n-1
    ADC = ADC + (abs(pitchingRate(idxADC))/maxPitchRate)*Ts;
end

ADC = ADC/(timeWindow(end) - timeWindow(1));

end