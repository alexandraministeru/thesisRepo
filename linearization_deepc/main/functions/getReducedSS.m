function LTIsys_reduced = getReducedSS(MBC,LTIsys,inputChannels,outputChannels)
inputChannelsList = MBC.DescCntrlInpt;
outputChannelsList = MBC.DescOutput;

nIn = length(inputChannels);
nOut = length(outputChannels);
nStates = length(LTIsys.A);

B_reduced = zeros(nStates,nIn);
C_reduced = zeros(nOut,nStates);
D_temp = zeros(length(LTIsys.C),nIn);
D_reduced = zeros(nOut,nIn);

for idx = 1:nIn
    % Find index of input channel
    id = find(ismember(inputChannelsList,inputChannels{idx}));
    B_reduced(:,idx) = MBC.AvgB(:,id);
    D_temp(:,idx) = MBC.AvgD(:,id);
end

for idx = 1:nOut
    % Find index of output channel
    id = find(ismember(outputChannelsList,outputChannels{idx}));
    C_reduced(idx,:) = MBC.AvgC(id,:);
    D_reduced(idx,:) = D_temp(id,:);
end

LTIsys_reduced = ss(MBC.AvgA,B_reduced,C_reduced,D_reduced);

end