function avg = getSSMean(data, ssWindowIdx, channels, channelName)
    dataIdx = find(ismember(channels,channelName));
    avg = mean(data(end-ssWindowIdx:end, dataIdx));
end