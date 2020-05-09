clf;close all;clear all

sig = [-1,1,-1,-1,1,-1,-1,-1,-1]

symbol_rate = 1*10^6
carrier_frequency = 16*10^6
IF_frequency = 4*10^6

srrc_2 = srrc_pulse(2, 10, 0.3);
srrc_16 = srrc_pulse(16,10, 0.3);

load('IIR_filter')

t_DAC_sig = conv(DAC(sig,16),srrc_2,'same')
t_DMA_sig = filter(IIR_filter,DAC(t_DAC_sig,4))


% modulatoin part
t_vec = [1:length(t_DMA_sig)]
carrier = cos(2*pi*(carrier_frequency-IF_frequency)/(32*10^6)*t_vec)
sig_with_carrier = t_DMA_sig .* carrier

t_DMA_sig = ADC(conv(sig_with_carrier,srrc_2,'same'),4)

% demodulatoin with carrier
t = [0:length(t_DMA_sig)-1]
carrier = exp(i*(2*pi*IF_frequency/(16*10^6)*t))
demod_sig = t_DMA_sig .* carrier

% demodulation
r_ADC_sig = ADC(conv(demod_sig,srrc_16,'same'),16)

r_sig = r_ADC_sig./max(abs(r_ADC_sig))

rcv_sig = (r_sig>0) + (-(r_sig<0))

subplot(2,1,1);stem(sig);title('BPSK Signal');grid on;
subplot(2,1,2);stem(rcv_sig);title('Recieved Signal');grid on;
figure()
plot(abs(fft(demod_sig)))

function DAC_sig = DAC(origin_sig,up_factor)
  DAC_sig = zeros(1,up_factor*length(origin_sig))
  DAC_sig([1:up_factor:length(DAC_sig)]) = origin_sig
end

function ADC_sig = ADC(origin_sig,down_factor)
  ADC_sig = zeros(1,length(origin_sig))
  ADC_sig = origin_sig(1:down_factor:end)
end

% SRRC Pulse
function [y,t] = srrc_pulse(T,A,a)
  t = [-A*T:A*T] + 10^(-8)
  if (a>0 && a<=1)
    num = cos((1+a)*pi*t/T) + T*sin((1-a)*pi*t/T)./(4*a*t)
    denom = 1-(4*a*t/T).^2
    y = 4*a/pi * num./denom
    y = y/max(y)
  else
    y = 1/T * sin(pi.*t./T) / (pi.*t./T)
    y = y/max(y)
  end
end

function y = add_awgn(x,SNR_DB)
  len = length(x)
  snr = 10^(SNR_DB/10)
  SYME = abs(x).^2 / len
  N0 = SYME / snr
  if isreal(x)
    n = sqrt(N0) * randn(1,len)
  else
    n = sqrt(N0/2) * (randn(1,len)+j*randn(1,len))
  y = x+n
  end
