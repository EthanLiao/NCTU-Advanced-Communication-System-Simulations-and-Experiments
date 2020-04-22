clf;clear all;
M = 16 % modulation order
k = log2(M) % modulation bits
n = 2000 %number of symbol
dataIn = randi([0 1],n,1) % generate a vector of data
dataInTupple = reshape(dataIn,n/k,k) % generate a group of tupple data
dataInSymbol = bi2de(dataInTupple) % transfer each tupple data into dec number
dataInSymbol = dataInSymbol.'
Gray_dataIn_recieve = QAMmod(dataInSymbol) % project data point to the QAM
SNR = 0:5:100
sym_err_array = zeros(1,length(SNR))
for idx=1:length(SNR)
  Gray_signal_recieve = add_awgn_noise(Gray_dataIn_recieve,SNR(idx))
  dataSymbolsOut = QAMdemod(Gray_signal_recieve)
  err_array = (dataSymbolsOut == dataInSymbol)
  symbol_error_rate = 1-sum(err_array)/length(err_array)
  sym_err_array(idx) = symbol_error_rate
end

% theoreitical symbol error rate
theoSymBer = 3/2*erfc(sqrt(0.1*(10.^(SNR/10))))
figure
plot(SNR,theoSymBer,'b.-','LineWidth',2)
hold on
plot(SNR,sym_err_array,'mx-','Linewidth',1)
axis = ([0 100 10^-5 1])
grid on
xlabel('ES/N0, dB')
ylabel('Symbol Error Rate')
legend('theory', 'simulation');
title('16-QAM Modulation : symbol error rate for theoreitical versus simulation');

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

function r_sig = QAMdemod(sig)
  r_sig = real(sig)                                % 接收訊號
  x = (0:15)                                       % 模擬16-QAM的每個訊號在複數平面上的投影
  std_QAM = QAMmod(x)
  for i = 1:length(sig)
    [argvalue, argmin] = min(abs(sig(i)-std_QAM))  % 將接收到的訊號與標準16-QAM的訊號相比，找距離最小的訊號
                                                   % 所代表的指標值
    r_sig(i) = x(argmin)                           % 以該指標值找出標準的16-QAM訊號
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
