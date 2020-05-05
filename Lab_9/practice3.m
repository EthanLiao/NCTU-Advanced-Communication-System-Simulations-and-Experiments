clf;close all;clear all

signal = [-1,1,-1,-1,1,-1,-1,-1,-1]

symbol_rate = 1*10^6
carrier_frequency = 8*10^6
IF_frequency = 4*10^6
srrc_2 = srrc_pulse(2, 10, 0.3);
srrc_16 = srrc_pulse(16,10, 0.3);

load('IIR_filter')
t_Digital_filtered_signal = conv(DAC(signal,16),srrc_2,'same')
t_DMA_filtered_signal = filter(IIR_filter,DAC(t_Digital_filtered_signal,2))


% modulatoin part
t_vec = [1:length(t_DMA_filtered_signal)]
carrier = cos(2*pi*(carrier_frequency-IF_frequency)/(32*10^6)*t_vec)

signal_with_carrier = t_DMA_filtered_signal .* carrier
t_DMA_filtered_signal = ADC(conv(signal_with_carrier,srrc_2,'same'),2)

% de modulatoin with carrier
t = [0:length(t_DMA_filtered_signal)-1]
carrier = exp(i*(2*pi*IF_frequency/(16*10^6)*t))
demodulation_signal = t_DMA_filtered_signal .* carrier

% demodulation
r_Digital_filtered_signal = ADC(conv(demodulation_signal,srrc_16,'same'),16)

r_sig = r_Digital_filtered_signal./max(abs(r_Digital_filtered_signal))
rcv_sig = (r_sig>0) + (-(r_sig<0))

subplot(2,1,1);stem(signal);title('BPSK Signal');grid on;
subplot(2,1,2);stem(rcv_sig);title('Recieved Signal');grid on;
figure()
plot(abs(fft(demodulation_signal)))

function DAC_sig = DAC(origin_signal,up_factor)
  DAC_sig = zeros(1,up_factor*length(origin_signal))
  DAC_sig([1:up_factor:length(DAC_sig)]) = origin_signal
end

function ADC_sig = ADC(origin_signal,down_factor)
  ADC_sig = zeros(1,length(origin_signal))
  ADC_sig = origin_signal(1:down_factor:end)
end

% SRRC Pulse
function [y,t] = srrc_pulse(T,A,a)
  t = [-A*T:A*T] + 10^(-8)
  if (a>0 && a<=1)
    num = cos((1+a)*pi*t/T) + T*sin((1-a)*pi*t/T)./(4*a*t)
    denom = 1-(4*a*t/T).^2
    y = 4*a/pi * num./denom
  else
    y = 1/T * sin(pi.*t./T) / (pi.*t./T)
  end
end
