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

% Generate 16-QAM sequence
M = 16                                      % modulation order
k = log2(M)                                 % modulation bits
n = 2000                                    % number of symbol
dataIn = randi([0 1],n,1)                   % generate a vector of data
dataInTupple = reshape(dataIn,n/k,k)        % generate a group of tupple data
dataInSymbol = bi2de(dataInTupple)          % transfer each tupple data into dec number
dataInSymbol = dataInSymbol.'
Gray_dataIn_recieve = QAMmod(dataInSymbol)  % project data point to the QAM

% SRRC Pulse
function [y,t] = srrc_pulse(T,A,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y = srrc_pulse(T, A, a)                                                       %
% OUTPUT                                                                        %
%      y: truncated SRRC pulse, with parameter T,                               %
%                 roll-off factor a, and duration 2*A*T                         %
%      t:   time axis of the truncated pulse                                    %
% INPUT                                                                         %
%      T:  Nyquist parameter or symbol period  (real number)                    %
%      A:  half duration of the pulse in symbol periods (positive INTEGER)      %
%      a:  roll-off factor (real number between 0 and 1)                        %
%                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  t = [-A*T:A*T] + 10^(-8)
  if (a>0 && a<=1)
    num = cos((1+a)*pi*t./T) + T*sin((1-a)*pi*t./T)./(4*a*t)
    denom = 1-(4*a*t./T).^2
    y = 4*a/pi * num./denom
  else
    y = 1/T * sin(pi.*t./T) / (pi.*t./T)
  end
end

% Gaussian Filter
function g_filter = Gfilter(BT, M, Tb)
  t = [-64:64];
  B = BT/Tb;
  C = sqrt(2*pi/log(2)) * B;
  g_filter = C*exp(-2*(pi^2)/log(2)*(BT/M)^2*(t.^2));
  g_filter = g_filter ./ sqrt(sum(g_filter.^2));
end
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
    n = sqrt(N0/2) * (randn(1,L)+i*randn(1,L))
  end
  y = x + n
end
% ----------- SNR ---------------------------
function snr = SNR(tr_sig,rcv_sig)
  if isreal(tr_sig)
    snr = 10*log10(mean(rcv_sig.^2)/mean((rcv_sig-tr_sig).^2));
  else
    snr = 10*log10(mean(abs(rcv_sig).^2)/mean(abs(rcv_sig-tr_sig).^2));
  end
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

% ------------Sampling Utility---------------------
function down_sig = ADC(signal,down)
  down_sig = zeros(1,length(signal))
  down_sig = signal(1:down:end)
end

function up_sig = DAC(signal,up)
  up_sig = zeros(1,up*length(signal))
  up_sig(1:up:end) = signal
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


function demod_sig = BPSK_demod(sig)
  idx = abs(sig)>0.3
  sig = sig(idx)
  pos_arr = sig>0
  neg_arr = sig<0
  neg_arr = neg_arr*(-1)
  demod_sig = pos_arr+neg_arr
end

function g_filter = Gfilter(BT, M, Tb)
  t = [-64:64];
  B = BT/Tb;
  C = sqrt(2*pi/log(2)) * B;
  g_filter = C*exp(-2*(pi^2)/log(2)*(BT/M)^2*(t.^2));
  g_filter = g_filter ./ sqrt(sum(g_filter.^2));
end
