M = 16 % modulation order
k = log2(M) % modulation bits
n = 10000 %number of symbol

dataIn = randi([0 1],n,1) % generate a vector of data
dataInTupple = reshape(dataIn,n/k,k) % generate a group of tupple data
dataInSymbol = bi2de(dataInTupple) % transfer each tupple data into dec number
Gray_dataIn_recieve = qammod(dataInSymbol,M) % project data point to the QAM



EbN0 = 10
snr = EbN0+10*log10(k)
Gray_signal_recieve = awgn(Gray_dataIn_recieve,snr,'measured') % generate noise around QAM data
noiseData = scatterplot(Gray_signal_recieve,1,0,'g.');
hold on;
scatterplot(Gray_dataIn_recieve,1,0,'k*',noiseData)
