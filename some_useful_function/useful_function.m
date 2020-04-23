% ----------- Some useful waves ---------------
% sinusoidal wave
fc = 3              % nature frequency of sinusoidal wave
N = 3               % number of sinusoidal wave
t = (0:0.01:N)/fc
signal_freq3 = cos(2*pi*fc*t)

% rectangular wave
rec_len = 10
mid = ceil(rec_len/2)
half = 3
rect = zeros(1,rec_len)
rect(mid-half:mid+half) = 1

% triangular wave
tri_wave = conv(rect,rect)
tri_wave_with_noise = add_awgn_noise(tri_wave,15)

% ------------- AWGN ------------------------------
function y = add_awgn_noise(x,SNR_DB)
  L = length(x)
  % calculate symbol energy
  SNR = 10^(SNR_DB/10) % SNR enery to linear scale
  SYME = sum(abs(x).^2) / L
  N0 = SYME / SNR      % Noise spectral Density
  if isreal(x)
    n = sqrt(N0) * randn(1,L)
  else
    n = (N0/2) * (randn(1,L)+i*randn(1,L))
  end
  y = x + n
end

% ----------- Padding Utility ---------------

function pad_array = pad(signal,zero_num)
  pad_array = zeros(1,(length(signal)+2*zero_num))
  pad_array((1+zero_num):(end-zero_num)) = signal(1:end)
end

function pad_array = pad_version2(signal,zero_num)
  pad_array = [zeros(1,zero_num) signal zeros(1,zero_num)]
end

function sh_signal = shift(signal,shamt)
  sh_signal = zeros(1,length(signal))
  sh_signal(1+shamt:end) = signal(1+shamt:end)
end

function delayed_signal = delay(sig,delay_num)
  delta = zeros(1,2*delay_num)
  delta(delay_num) = 1
  delayed_signal = conv(sig,delta)
end

function delayed_signal = delay_version2(sig,delay_num)
  delayed_signal = [zeros(1,delay_num) sig]
end

% -------------Modulation--------------------------
function c_sig = PAMmod(sig)
  for i=1:length(sig)
    if sig(i) == 1
      c_sig(i) = -7
      continue
    elseif sig(i) == 2
      c_sig(i) = -5
      continue
    elseif sig(i) == 3
      c_sig(i) = -3
      continue
    elseif sig(i) == 4
      c_sig(i) = -1
      continue
    elseif sig(i) == 5
      c_sig(i) = 1
      continue
    elseif sig(i) == 6
      c_sig(i) = 3
      continue
    elseif sig(i) == 7
      c_sig(i) = 5
      continue
    elseif sig(i) == 8
      c_sig(i) = 7
      continue
    end
  end
end

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
    elseif sig(i)==complex(9) % 10
      c_sig(i) = 3+j
      continue
    elseif sig(i)==complex(10) % 10
      c_sig(i) = 3-3*j
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
