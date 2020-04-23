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

function c_sig = QAMmod(sig)
  sig = complex(sig)
  for i=1:length(sig)
    if sig(i)== complex(0)     % 0
      c_sig(i) = -3+3*j
      continue
    elseif sig(i)==complex(1) % 1
      c_sig(i) = -3+j
      continue
    elseif sig(i)==complex(2) % 2
      c_sig(i) = -3-3*j
      continue
    elseif sig(i)==complex(3) % 3
      c_sig(i) = -3-j
      continue
    elseif sig(i)==complex(4) % 4
      c_sig(i) = -1+3*j
      continue
    elseif sig(i)==complex(5) % 5
      c_sig(i) = -1+j
      continue
    elseif sig(i)==complex(6) % 6
      c_sig(i) = -1-3*j
      continue
    elseif sig(i)==complex(7) % 7
      c_sig(i) = -1-j
      continue
    elseif sig(i)==complex(8) % 8
      c_sig(i) = 3+3*j
      continue
    elseif sig(i)==complex(10) % 10
      c_sig(i) = 3-3*j
      continue
    elseif sig(i)==complex(9) % 10
      c_sig(i) = 3+j
      continue
    elseif sig(i)==complex(11) % 11
      c_sig(i) = 3-j
      continue
    elseif sig(i)==complex(12) % 12
      c_sig(i) = 1+3*j
      continue
    elseif sig(i)==complex(13) % 13
      c_sig(i) = 1+j
      continue
    elseif sig(i)==complex(14) % 14
      c_sig(i) = 1-3*j
      continue
    elseif sig(i)==complex(15) % 15
      c_sig(i) = 1-j
      continue
    end
  end
end
