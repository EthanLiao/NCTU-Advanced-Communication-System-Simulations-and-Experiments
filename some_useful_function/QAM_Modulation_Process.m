clf;clear all;

M = 16
k = log2(M)
N = 2000
SNR_db = 20
x = (0:15)
std_QAM = QAMmod(x)
bit_sig = randi([0,1],N,1)
bin_sig = reshape(bit_sig,N/k,k)
symbol_sig = bi2de(bin_sig).'
QAM_sig = QAMmod(symbol_sig)
awgn_sig = add_awgn_noise(QAM_sig,SNR_db)
demod_sig = QAMdmod(awgn_sig)

symERR = 1-sum(demod_sig==symbol_sig)/length(symbol_sig)

noise_data = scatterplot(awgn_sig,1,0,'g.')
hold on;
scatterplot(std_QAM,1,0,'k*',noise_data)

function orig_sig = QAMdmod(sig)
  x = (0:15)
  std_QAM = QAMmod(x)
  for i = 1:length(sig)
    [minvalue,minarg] = min(abs(sig(i)-std_QAM))
    orig_sig(i) = x(minarg)
  end
end




function QAM_sig = QAMmod(sig)
  for i=1:length(sig)
    if sig(i) == 0
      QAM_sig(i) = -3+3*j
      continue
    elseif sig(i) == 1
      QAM_sig(i) = -3+1*j
      continue
    elseif sig(i) == 2
      QAM_sig(i) = -3-3*j
      continue
    elseif sig(i) == 3
      QAM_sig(i) = -3-1*j
      continue
    elseif sig(i) == 4
      QAM_sig(i) = -1+3*j
      continue
    elseif sig(i) == 5
      QAM_sig(i) = -1+1*j
      continue
    elseif sig(i) == 6
      QAM_sig(i) = -1-3*j
      continue
    elseif sig(i) == 7
      QAM_sig(i) = -1-1*j
      continue
    elseif sig(i) == 8
      QAM_sig(i) = 3+3*j
      continue
    elseif sig(i) == 9
      QAM_sig(i) = 3+1*j
      continue
    elseif sig(i) == 10
      QAM_sig(i) = 3-3*j
      continue
    elseif sig(i) == 11
      QAM_sig(i) = 3-1*j
      continue
    elseif sig(i) == 12
      QAM_sig(i) = 1+3*j
      continue
    elseif sig(i) == 13
      QAM_sig(i) = 1+1*j
      continue
    elseif sig(i) == 14
      QAM_sig(i) = 1-3*j
      continue
    elseif sig(i) == 15
      QAM_sig(i) = 1-1*j
      continue
    end
  end
end

function y = add_awgn_noise(x,SNR_db)
  rng('default')
  L = length(x)
  snr = 10^(SNR_db/10)
  SYME = sum(abs(x).^2)/L
  N0 = SYME/snr
  if isreal(x)
    n = sqrt(N0) * randn(1,L)
  else
    n = (N0/2) * (randn(1,L)+i*randn(1,L))
  end
  y = x+n
end
