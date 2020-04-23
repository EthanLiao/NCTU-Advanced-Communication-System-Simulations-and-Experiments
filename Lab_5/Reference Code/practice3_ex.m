clf
M = 16 % modulation order
k = log2(M) % modulation bits
n = 3000 %number of bit to process
data = 4
q_data = qammod(data,M)
x = (0:15)' % generate a modulation symbol
y = qammod(x,M,'gray')
[argvalue, argmin] = min(abs(y.'-q_data))
ans = x(argmin)



M = 16 % modulation order
k = log2(M) % modulation bits
n = 100000 %number of symbol
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
% convert recived signal to binary
dataSymbolsOut = qamdemod(Gray_signal_recieve,M,'bin')
dataOutMatrix = de2bi(dataSymbolsOut,k)
dataOut = dataOutMatrix(:)
[numErrors,ber] = biterr(dataIn,dataOut);
[number,ratio] = symerr(dataInSymbol,dataSymbolsOut)


function num = GrayToBinary(num)
   mask = bitshift(num, -1);
   while mask > 0
      num = bitxor(num, mask);
      mask = bitshift(mask, -1);
   end
end

function num = BinaryToGray(num)
  if 0000     % 0
    num = 0000
  elseif 0001 % 1
    num = 0001
  elseif 0010 % 2
    num = 0011
  elseif 0011 % 3
    num = 0010
  elseif 0100 % 4
    num = 0110
  elseif 0101 % 5
    num = 0111
  elseif 0110 % 6
    num = 0101
  elseif 0111 % 7
    num = 0100
  elseif 1000 % 8
    num = 1100
  elseif 1001 % 9
    num = 1101
  elseif 1010 % 10
    num = 1111
  elseif 1011 % 11
    num = 1110
  elseif 1100 % 12
    num = 1010
  elseif 1101 % 13
    num = 1011
  elseif 1110 % 14
    num = 1001
  elseif 1111 % 15
    num = 1000
  end
end
