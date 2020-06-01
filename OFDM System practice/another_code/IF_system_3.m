clf;close all;clear all
N = 20;
sig = randi([0,1],1,N);
sig((sig==0)) = -1;

symbol_rate = 1*10^6;
carrier_frequency = 32*10^6;
IF_frequency = 2*10^6;
freq_DAC = 16*10^6;
freq_DMA = 128*10^6;
f_DAC = 16;
f_DMA = 8;

srrc_16 = srrc_pulse(16, 5, 0.3);
srrc_8 = srrc_pulse(2,5,0.3);

srrc_16_delay = (length(srrc_16)-1)/2;
srrc_8_delay = (length(srrc_8)-1)/2;

load('./filter/IIR_filter')

% modulatoin part
t_DAC_sig = conv(DAC(sig,f_DAC),srrc_16);
t_DAC_sig = t_DAC_sig(srrc_16_delay+1:end-srrc_16_delay);

t_DMA_ana_sig = conv(DAC(t_DAC_sig,f_DMA),srrc_8);
t_DMA_ana_sig = t_DMA_ana_sig(srrc_8_delay+1:end-srrc_8_delay);

t_DMA_sig = filter(IIR_filter,t_DMA_ana_sig);

% modulatoin with carrier
t = [0:length(t_DMA_sig)-1];
carrier = sqrt(2)*exp(1j*2*pi*carrier_frequency*t);
sig_with_carrier = real(t_DMA_sig .* carrier);

% IF Band
t_vec = [1:length(sig_with_carrier)];
carrier = cos(2*pi*(carrier_frequency-IF_frequency)/freq_DMA*t_vec);

IF_sig = sig_with_carrier .* carrier;
r_DMA_f_sig = conv(IF_sig,srrc_8);
r_DMA_f_sig = r_DMA_f_sig(srrc_8_delay+1:end-srrc_8_delay);
r_DMA_sig = ADC(r_DMA_f_sig,f_DMA);

% demodulatoin with carrier
t = [0:length(r_DMA_sig)-1];
carrier = exp(-1j*(2*pi*IF_frequency/freq_DAC*t+pi));
demod_sig = real(r_DMA_sig .* carrier);

% demodulation
r_ADC_f_sig = conv(demod_sig,srrc_16);
r_ADC_f_sig = r_ADC_f_sig(srrc_16_delay+1:end-srrc_16_delay);
r_ADC_sig = ADC(r_ADC_f_sig,f_DAC);

r_sig = r_ADC_sig./max(abs(r_ADC_sig));

% signal detection
% rcv_sig = (r_sig>0) + (-(r_sig<0));
snr = SNR(sig, r_sig)

subplot(2,1,1);stem(sig);title('BPSK Signal');grid on;
subplot(2,1,2);stem(r_sig);title('Recieved Signal');grid on;
figure()
plot(abs(fft(demod_sig)));

function DAC_sig = DAC(origin_sig,up_factor)
  DAC_sig = zeros(1,up_factor*length(origin_sig));
  DAC_sig([1:up_factor:length(DAC_sig)]) = origin_sig;
end

function ADC_sig = ADC(origin_sig,down_factor)
  ADC_sig = zeros(1,length(origin_sig));
  ADC_sig = origin_sig(1:down_factor:end);
end

% SRRC Pulse
function [y,t] = srrc_pulse(T,A,a)
  t = [-A*T:A*T] + 10^(-8);
  if (a>0 && a<=1)
    num = cos((1+a)*pi*t/T) + T*sin((1-a)*pi*t/T)./(4*a*t);
    denom = 1-(4*a*t/T).^2;
    y = 4*a/pi * num./denom;
  else
    y = 1/T * sin(pi.*t./T) ./ (pi*t./T);
  end
end

function snr = SNR(tr_sig,rcv_sig)
  if isreal(tr_sig)
    snr = mean(rcv_sig.^2)/mean((rcv_sig-tr_sig).^2);
  else
    snr = mean(abs(rcv_sig).^2)/mean(abs(rcv_sig-tr_sig).^2);
  end
end
