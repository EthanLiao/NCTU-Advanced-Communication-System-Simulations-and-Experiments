clf;close all;clear all
N = 20;
sig = randi([0 1],1,N);
sig((sig==0)) = -1;

symbol_rate = 1e6;
carrier_frequency = 16e6;
IF_frequency = 4e6;
freq_DAC = 16e6;
freq_DMA = 64e6;
f_DAC = freq_DAC / symbol_rate;
f_DMA = freq_DMA / freq_DAC;

SNR_DB = 20;
f_noise = 0.9;
srrc_16 = srrc_pulse(16, 5, 0.3);
srrc_4 = srrc_pulse(4,5,0.3);

srrc_16_delay = (length(srrc_16)-1)/2;
srrc_4_delay = (length(srrc_4)-1)/2;

load('./filter/IIR_64_filter')
load('./filter/clr_filter')

group_delay_IIR = 9;
group_delay_clr = 1;
% modulatoin part

t_DAC_sig = conv(DAC(sig,f_DAC),srrc_16);
t_DAC_sig = t_DAC_sig(srrc_16_delay+1:end-srrc_16_delay);

t_DMA_ana_sig = conv(DAC(t_DAC_sig,f_DMA),srrc_4);
t_DMA_ana_sig = t_DMA_ana_sig(srrc_4_delay+1:end-srrc_4_delay);

t_DMA_sig = filter(IIR_64_filter,t_DMA_ana_sig);
t_DMA_sig = t_DMA_sig(group_delay_IIR:end);

% modulatoin with carrier
t = [0:length(t_DMA_sig)-1];
carrier = sqrt(2)*exp(1j*2*pi*carrier_frequency/freq_DMA*t);
sig_with_carrier = real(t_DMA_sig .* carrier);

% awgn_sig = add_awgn_noise(sig_with_carrier,SNR_DB);
awgn_sig = awgn(sig_with_carrier,SNR_DB);
t = [0:length(awgn_sig)-1];
sig_noise = cos(2*pi*f_noise*t);

% Image rejection Filter
dirty_sig = awgn_sig+sig_noise;
clear_sig = filter(clr_filter,dirty_sig);
clear_sig = clear_sig(group_delay_clr:end);
% clear_sig = awgn_sig;
% dirty_sig = awgn_sig;
% IF Band
t_vec = [0:length(clear_sig)-1];
carrier = cos(2*pi*(carrier_frequency-IF_frequency)/freq_DMA*t_vec);
IF_sig = clear_sig .* carrier;

t_vec = [0:length(dirty_sig)-1];
carrier = cos(2*pi*(carrier_frequency-IF_frequency)/freq_DMA*t_vec);
IF_dirty_sig = dirty_sig.*carrier;

r_DMA_f_sig = conv(IF_sig,srrc_4);
r_DMA_f_sig = r_DMA_f_sig(srrc_4_delay+1:end-srrc_4_delay);
r_DMA_sig = ADC(r_DMA_f_sig,f_DMA);

% demodulatoin with carrier
t = [0:length(r_DMA_sig)-1];
carrier = exp(-1j*2*pi*IF_frequency/freq_DAC*t);
demod_sig = real(r_DMA_sig .* carrier);

% demodulation
r_ADC_f_sig = conv(demod_sig,srrc_16);
r_ADC_f_sig = r_ADC_f_sig(srrc_16_delay+1:end-srrc_16_delay);
r_ADC_sig = ADC(r_ADC_f_sig,f_DAC);

r_sig = r_ADC_sig./max(r_ADC_sig);
rcv_sig = (r_sig>0) + (-(r_sig<0));

subplot(3,2,1);stem(sig);title('Transmitting Signal');grid on;
subplot(3,2,2);plot(abs(fft(sig)));title('Transmitting Signal fdomain');grid on;
subplot(3,2,3);stem(IF_dirty_sig);title('IF Band Interference Signal');grid on;
subplot(3,2,4);plot(abs(fft(IF_dirty_sig)));title('IF Band Interference Signal fdomain');grid on;
subplot(3,2,5);stem(rcv_sig);title('Recieved Signal');grid on;
subplot(3,2,6);plot(abs(fft(rcv_sig)));title('Recieved Signal fdomain');grid on;

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

function y = add_awgn_noise(x,SNR_DB)
  snr = 10^(SNR_DB/10);
  L = length(x);
  SYME = sum(abs(x).^2) / L;
  N0 = SYME / snr;
  if isreal(x)
    n = sqrt(N0)*randn(1,L);
  else
    n = sqrt(N0/2)*(randn(1,L)+j*randn(1,L));
  end
  y = x+n;
end
