clf;clear all;
M = 16                                    % modulation order
k = log2(M)                                 % modulation bits
n = 2000                                    % number of symbol
dataIn = randi([0 1],n,1)                   % generate a vector of data
dataInTupple = reshape(dataIn,n/k,k)        % generate a group of tupple data
dataInSymbol = bi2de(dataInTupple)          % transfer each tupple data into dec number
dataInSymbol = dataInSymbol.'
Gray_dataIn_recieve = QAMmod(dataInSymbol)  % project data point to the QAM
SNR = 10
Gray_signal_recieve = add_awgn_noise(Gray_dataIn_recieve,SNR)
dataSymbolsOut = QAMdemod(Gray_signal_recieve)
err_array = dataSymbolsOut == dataInSymbol
symbol_error_rate = 1-sum(err_array)/length(err_array)

noiseData = scatterplot(Gray_signal_recieve,1,0,'g.');
hold on;
scatterplot(Gray_dataIn_recieve,1,0,'k*',noiseData)
title(['symbol error rate : ' num2str(symbol_error_rate)] )

function c_sig = QAMmod(sig)
  sig = complex(sig)
  for i=1:length(sig)
    if sig(i)== complex(0)    % 0
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

function r_sig = QAMdemod(sig)
  r_sig = real(sig)
  x = (0:15)
  std_QAM = QAMmod(x)
  for i = 1:length(sig)
    [argvalue, argmin] = min(abs(sig(i)-std_QAM))
    r_sig(i) = x(argmin)
  end
end

function y = add_awgn_noise(x,SNR_dB)
 %y=awgn_noise(x,SNR) adds AWGN noise vector to signal 'x' to generate a
 %resulting signal vector y of specified SNR in dB
 rng('default');%set the random generator seed to default (for comparison only)
 L=length(x);
 SNR = 10^(SNR_dB/10); %SNR to linear scale
 Esym=sum(abs(x).^2)/(L); %Calculate actual symbol energy
 N0=Esym/SNR; %Find the noise spectral density
 if(isreal(x)),
 noiseSigma = sqrt(N0);%Standard deviation for AWGN Noise when x is real
 n = noiseSigma*randn(1,L);%computed noise
 else
 noiseSigma=sqrt(N0/2);%Standard deviation for AWGN Noise when x is complex
 n = noiseSigma*(randn(1,L)+1i*randn(1,L));%computed noise
 end
 y = x + n; %received signal
end
