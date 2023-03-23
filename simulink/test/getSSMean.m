function avg = getSSMean(data, ssWindowIdx, dataIdx)
    avg = mean(data(end-ssWindowIdx:end, dataIdx));
end