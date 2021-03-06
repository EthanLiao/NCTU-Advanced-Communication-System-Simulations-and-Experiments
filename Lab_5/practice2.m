M = 16                                % modulation order
k = log2(M)                           % modulation bits
n = 10000                             %number of symbol

dataIn = randi([0 1],n,1)             % generate a vector of data
dataInTupple = reshape(dataIn,n/k,k)  % generate a group of tupple data
dataInSymbol = bi2de(dataInTupple)    % transfer each tupple data into dec number
Gray_dataIn_recieve = QAMmod(dataInSymbol)
Gray_signal_recieve = add_awgn_noise(Gray_dataIn_recieve,10)
noiseData = scatterplot(Gray_signal_recieve,1,0,'g.');
hold on;
scatterplot(Gray_dataIn_recieve,1,0,'k*',noiseData)

function y = add_awgn_noise(x,snr_db)
  L = length(x)
  SNR = 10^(snr_db/10)
  SYME = sum(abs(x.^2))/L
  N0 = SYME/SNR
  if isreal(x)
    n = sqrt(N0)*randn(1,L)
  else
    n = (N0/2) * (randn(1,L)+i*randn(1,L))
  end
  y = x + n
end

function data_modulated = QAM_mod(dec_data)
  for i=1:length(dec_data)
    if dec_data(i) == 0
      data_modulated(i) = -3+3*j
      continue
    elseif dec_data(i) == 1
      data_modulated(i) = -3+1*j
      continue
    elseif dec_data(i) == 2
      data_modulated(i) = -3-1*j
      continue
    elseif dec_data(i) == 3
      data_modulated(i) = -3-3*j
      continue
    elseif dec_data(i) == 4
      data_modulated(i) = -1+3*j
      continue
    elseif dec_data(i) == 5
      data_modulated(i) = -1+1*j
      continue
    elseif dec_data(i) == 6
      data_modulated(i) = -1-3*j
      continue
    elseif dec_data(i) == 7
      data_modulated(i) = -1-1*j
      continue
    elseif dec_data(i) == 8
      data_modulated(i) = 3+3*j
      continue
    elseif dec_data(i) == 9
      data_modulated(i) = 3+1*j
      continue
    elseif dec_data(i) == 10
      data_modulated(i) = 3-3*j
      continue
    elseif dec_data(i) == 11
      data_modulated(i) = 3-1*j
      continue
    elseif dec_data(i) == 12
      data_modulated(i) = 1+3*j
      continue
    elseif dec_data(i) == 13
      data_modulated(i) = 1+1*j
      continue
    elseif dec_data(i) == 14
      data_modulated(i) = 1-3*j
      continue
    else dec_data(i) == 15
      data_modulated(i) = 1-1*j
      continue
    end
  end
end
