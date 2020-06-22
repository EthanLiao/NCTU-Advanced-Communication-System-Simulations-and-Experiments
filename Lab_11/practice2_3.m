clf;clear all;close all;

signal = [1 zeros(1,127)];

% Causal System
ADD_AWGN = false;
SNR_DB = 0;
channel_gain_fast = 1 ; channel_gain_slow = -0.5;
ISI_signal = real(ISI_system(signal, channel_gain_fast, channel_gain_slow, ADD_AWGN, SNR_DB));
L = length(ISI_signal);
h = 0.5.^[0:L-1];
equalized_signal = conv(ISI_signal,h);
equalized_signal = equalized_signal / max(equalized_signal);
% Non - casaul System
ADD_AWGN = false;
SNR_DB = 0;
channel_gain_fast = 0.5 ; channel_gain_slow = -1;
ISI_signal = real(ISI_system(signal, channel_gain_fast, channel_gain_slow, ADD_AWGN, SNR_DB));
w = 2*(-2)*2.^[-L:-1];
causal_equalized_signal  = conv(ISI_signal,w);
causal_equalized_signal = causal_equalized_signal / max(causal_equalized_signal);

% ADD AWGN : 5 DB in non casaul system
ADD_AWGN = true;
SNR_DB = -10;
channel_gain_fast = 0.5 ; channel_gain_slow = -1;
ISI_signal = real(ISI_system(signal, channel_gain_fast, channel_gain_slow, ADD_AWGN, SNR_DB));
w = (-2)*2.^[-L:-1];
awgn_causal_equalized_signal  = conv(ISI_signal,w);
awgn_causal_equalized_signal = awgn_causal_equalized_signal/ max(awgn_causal_equalized_signal);

subplot(5,1,1);stem(signal);title('Transmitted Signal');grid on;
subplot(5,1,2);stem(ISI_signal);title('Recieced signal');grid on;
subplot(5,1,3);stem(equalized_signal);title('Recover Causal');grid on;
subplot(5,1,4);stem(causal_equalized_signal);title('Causal Recover Causal');grid on;
subplot(5,1,5);stem(awgn_causal_equalized_signal);title('AWGN Causal Recover Causal');grid on;

function isi_sig = ISI_system(signal, channel_gain_fast, channel_gain_slow, ADD_AWGN, SNR_DB)
  t_sig = trans_branch(signal);
  multipath_effect = [channel_gain_fast zeros(1, 63) channel_gain_slow];
  c_sig = conv(t_sig,multipath_effect);
  if ADD_AWGN
    c_sig = add_awgn_noise(c_sig, SNR_DB);
  end
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

function y = add_awgn_noise(x,SNR_DB)
  L = length(x);
  % calculate symbol energy
  SNR = 10^(SNR_DB/10); % SNR enery to linear scale
  SYME = sum(abs(x).^2) / L;
  N0 = SYME / SNR;      % Noise spectral Density
  if isreal(x)
    n = sqrt(N0) * randn(1,L);
  else
    n = sqrt(N0/2) * (randn(1,L)+i*randn(1,L));
  end
  y = x + n;
end
