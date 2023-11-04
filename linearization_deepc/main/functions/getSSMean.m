function avg = getSSMean(data, ssWindowIdx, channels, channelName)
%getSSMean(data, ssWindowIdx, channels, channelName) Function used to get
%the mean of the final data points on a specified channel. Given that the
%simulated model reached steady state within the simulation time, the
%function returns the steady-state value of the channel.
% 
% Input arguments
%----------------
% data          : data matrix read from a FAST .out file.
% ssWindowIdx   : width of the average window, the mean is calculated over
%                 the ssWindowIdx data points (so equal to time window (in
%                 seconds) multiplied by the sampling period (in seconds)).
% channels      : channels list read from a FAST .out file.
% channelName   : the name of the channel for which the mean will be
%                 computed, must be in the channels list above.
%
% Output arguments:
%------------------
% avg           : scalar representing the mean of the channel over the
%                 last ssWindowIdx data points.
%==========================================================================

dataIdx = find(ismember(channels,channelName));
avg = mean(data(end-ssWindowIdx:end, dataIdx)); %#ok<FNDSB>

end