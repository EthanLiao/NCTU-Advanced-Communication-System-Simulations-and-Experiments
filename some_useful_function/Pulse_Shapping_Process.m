clf;clear all;close all
N = 100
mid = ceil(N/2)
half = 20
UP_FACTOR = 4
DOWN_FACTOR = 4
SNR_DB = 15
rect = zeros(1,N)
rect(mid-half:mid+half) = 1
tri_pulse = conv(rect,rect)
srrc_4 = srrc_pulse(4,100,0.3)
srrc_16 = srrc_pulse(16,100,0.3)
% signal = tri_pulse

signal = [-1,1,1,-1,1,-1,1]

load('filter/pulse_shapping_filter')
%%%%%%%%% NON-Pratical SRRC Pulse Shapping %%%%%%%%%%%%

trans_sig = conv(DAC(signal,UP_FACTOR),srrc_4,'same')
rcv_sig = ADC(conv(trans_sig,srrc_4,'same'),DOWN_FACTOR)./100

stem(signal)
hold on;stem(rcv_sig)
title('Non Pratical SRRC Pulse Shapping')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Practical SRRC Pulse Shapping %%%%%%%%%%%%

trans_sig = conv(DAC(conv(DAC(signal,UP_FACTOR),srrc_4,'same'),UP_FACTOR),srrc_16)
awgn_sig = add_awgn_noise(trans_sig,SNR_DB)
rcv_sig = ADC(conv(ADC(conv(awgn_sig,srrc_16,'same'),DOWN_FACTOR),srrc_4,'same'),DOWN_FACTOR)

% calculate symbol error rate
rcv_demod_sig = BPSK_demod(rcv_sig./max(rcv_sig))
sym_ERR = 1-(sum(signal==rcv_demod_sig)./length(signal))

figure()
stem(delay(signal,100))
hold on;stem(rcv_sig./max(rcv_sig))
title(['Practical SRRC Pulse Shapping : ' num2str(sym_ERR)] )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Practical General Filter Pulse Shapping %%%%%%%%%%%%
trans_sig = filter(pulse_shapping_filter,DAC(conv(DAC(signal,UP_FACTOR),srrc_4,'same'),UP_FACTOR))
% awgn_sig = add_awgn_noise(trans_sig,SNR_DB)
rcv_sig = ADC(conv(ADC(filter(pulse_shapping_filter,trans_sig),DOWN_FACTOR),srrc_4,'same'),DOWN_FACTOR)

figure()
stem(delay(signal,0))
hold on;stem(rcv_sig./max(abs(rcv_sig)))
title('Practical General Filter Pulse Shapping , SER : ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function up_sig = DAC(sig,up_factor)
  up_sig = zeros(1,length(sig)*up_factor)
  up_sig(1:up_factor:end) = sig
end

function down_sig = ADC(sig,down_factor)
  down_sig = zeros(1,length(sig))
  down_sig = sig(1:down_factor:end)
end

function srrc = srrc_pulse(T,A,a)
  t = [-A*T:A*T]+10^-8
  if a>0 && a<=1
    num = cos((1+a)*t/T)+T*sin((1-a)*t/T)./(4*a*t/T)
    denom = 1-(4*a*t/T).^2
    srrc = 4*a/pi *num ./denom
  else
    srrc = 1/T * sin(pi*t/T) ./ (t/T)
  end
end

function y = add_awgn_noise(x,SNR_DB)
  L = length(x)
  snr = 10^(SNR_DB/10)
  SYME = sum(abs(x).^2)/L
  N0 = SYME/L
  if isreal(x)
    n = N0 * randn(1,L)
  else
    n = (N0/2) * (randn(1,L)+i*randn(1,L))
  end
  y = x+n
end

function demod_sig = BPSK_demod(sig)
  idx = abs(sig)>0.3
  sig = sig(idx)
  pos_arr = sig>0
  neg_arr = sig<0
  neg_arr = neg_arr*(-1)
  demod_sig = pos_arr+neg_arr
end

function pad_sig = padding(sig,padamt)
  pad_sig = [zeros(1,padamt) sig zeros(1,padamt)]
end

function delay_sig = delay(sig,delay)
  delay_sig = [zeros(1,delay) sig]
end
