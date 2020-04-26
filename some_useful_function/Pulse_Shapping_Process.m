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

signal = [-1,1,1,1]

%%%%%%%%% NON-Pratical SRRC Pulse Shapping %%%%%%%%%%%%

trans_sig = conv(ADC(signal,UP_FACTOR),srrc_4,'same')
rcv_sig = DAC(conv(trans_sig,srrc_4,'same'),DOWN_FACTOR)./100

subplot(3,1,1);stem(signal)
subplot(3,1,2);plot(trans_sig)
subplot(3,1,3);stem(rcv_sig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Practical SRRC Pulse Shapping %%%%%%%%%%%%

trans_sig = conv(ADC(conv(ADC(signal,UP_FACTOR),srrc_4,'same'),UP_FACTOR),srrc_16)
awgn_sig = add_awgn_noise(trans_sig,SNR_DB)
rcv_sig = DAC(conv(DAC(conv(awgn_sig,srrc_16,'same'),DOWN_FACTOR),srrc_4,'same'),DOWN_FACTOR)

figure()
stem(delay(signal,100))
hold on;stem(rcv_sig./max(rcv_sig))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Practical General Filter Pulse Shapping %%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function up_sig = ADC(sig,up_factor)
  up_sig = zeros(1,length(sig)*up_factor)
  up_sig(1:up_factor:end) = sig
end

function down_sig = DAC(sig,down_factor)
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

function delay_sig = delay(sig,delay)
  delay_sig = [zeros(1,delay) sig]
end
