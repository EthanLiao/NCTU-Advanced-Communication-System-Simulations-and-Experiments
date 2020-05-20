clf;clear all;close all;

signal = [1 zeros(1,127)];

% Causal System
channel_gain_fast = 0.3 ; channel_gain_med = 1 ;  channel_gain_slow = 0.3;
ISI_signal = real(ISI_system(signal, channel_gain_fast, channel_gain_med, channel_gain_slow));
L = length(ISI_signal);
h = 8*(-1/3)*((-1/3).^[0:L-1]);
equalized_signal = conv(ISI_signal,h);
L = length(equalized_signal);
w = 80*((-3).^[-L:-1]);
equalized_signal = conv(equalized_signal,w);

subplot(3,1,1);stem(signal);title('Transmitted Signal');grid on;
subplot(3,1,2);stem(ISI_signal);title('Recieced signal');grid on;
subplot(3,1,3);stem(equalized_signal);title('Recover Causal');grid on;

function isi_sig = ISI_system(signal, channel_gain_fast, channel_gain_med, channel_gain_slow)
  t_sig = trans_branch(signal);
  multipath = [channel_gain_fast zeros(1, 63) channel_gain_med zeros(1, 63) channel_gain_slow ];
  c_sig = conv(t_sig,multipath);
  isi_sig = recieve_branch(c_sig);
end

function trans_sig = trans_branch(sig)
  % modulatoin part
  fc = 16*10^6;
  fs = 64*10^6;
  f_DAC = 16;
  f_DMA = 4;
  srrc_4 =srrc_pulse(4, 5, 1);
  srrc_16 =srrc_pulse(16, 5, 1);
  t_DAC_sig = conv(DAC(sig, f_DAC), srrc_16, 'same');
  t_DMA_sig = conv(DAC(t_DAC_sig, f_DMA), srrc_4, 'same');
  t = [0:length(t_DMA_sig)-1];
  trans_sig = real(t_DMA_sig .* exp(j*2*pi*fc/fs*t));
end

function rcv_sig = recieve_branch(demod_sig)
  fc = 16*10^6;
  fs = 64*10^6;
  f_DAC = 16;
  f_DMA = 4;
  srrc_4 =srrc_pulse(4, 5, 1);
  srrc_16 =srrc_pulse(16, 5, 1);

  t = [0:length(demod_sig)-1];
  demod_sig = demod_sig .* exp(-j*2*pi*fc/fs*t);

  f_sig = conv(demod_sig, srrc_4, 'same');
  r_DMA_sig = ADC(f_sig,f_DMA);

  r_DAC_sig = conv(r_DMA_sig,srrc_16, 'same');
  rcv_sig = ADC(r_DAC_sig,f_DAC);
end

function DAC_sig = DAC(origin_signal,up_factor)
  DAC_sig = zeros(1,up_factor*length(origin_signal));
  DAC_sig([1:up_factor:length(DAC_sig)]) = origin_signal;
end

function ADC_sig = ADC(origin_signal,down_factor)
  ADC_sig = zeros(1,length(origin_signal));
  ADC_sig = origin_signal([1:down_factor:end]);
end

function [phi, t] = srrc_pulse(T, A, a)
t = [-A*T:A*T] + 10^(-8); % in order to avoid division by zero problems at t=0.
  if (a>0 && a<=1)
     num = cos((1+a)*pi*t/T) + sin((1-a)*pi*t/T) ./ (4*a*t/T);
     denom = 1-(4*a*t./T).^2;
     phi = 4*a/(pi*sqrt(T)) * num ./ denom;
  elseif (a==0)
     phi = 1/(sqrt(T)) * sin(pi*t/T)./(pi*t/T);
  end
end
