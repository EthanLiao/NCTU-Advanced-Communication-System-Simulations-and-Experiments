clf;clear all;close all;

% Causal System
imp = [1 zeros(1,49)];
zero = 1;
pole = [1, -1.2, -0.25, 0.3];
h = filter(zero, pole, imp);
stem(h);

tap = zeros(1,19);
N = 10;
% signal = [1 tap 1 tap 1 tap];
signal = [];
for i = 1:3
  signal = [signal 1 tap];
end

SNR_DB = 20;
multipath = [1 tap -1.2 tap -0.25 tap 0.3];

ADD_AWGN = false;
ISI_signal = ISI_system(signal, multipath, ADD_AWGN, SNR_DB);

L = length(ISI_signal);
h = (-20/56)*((1/2).^[0:L-1]);
equalized_signal = conv(ISI_signal,h);

L = length(equalized_signal);
w = 20*36/(119*5)*(6/5)*((6/5).^[-L:-1]);
equalized_signal = conv(equalized_signal,w);

L = length(equalized_signal);
w = 20/(68*2)*((-1/2).^[0:L-1]);
equalized_signal = conv(equalized_signal,w);


% figure()
% stem(equalized_signal);
group_delay = 126;
equalized_signal = equalized_signal(group_delay:end);


equalized_signal = equalized_signal(1:60);
equalized_signal = equalized_signal / max(equalized_signal);

snr_non_AWGN = SNR(signal,equalized_signal)

figure();
subplot(3,1,1);stem(signal);title('Transmitted Signal');grid on;
subplot(3,1,2);stem(ISI_signal);title('Recieced signal');grid on;
subplot(3,1,3);stem(equalized_signal);title('Recover Causal');grid on;


% add awgn
ADD_AWGN = true;
ISI_signal = real(ISI_system(signal, multipath, ADD_AWGN, SNR_DB));

L = length(ISI_signal);
h = (-20/56)*((1/2).^[0:L-1]);
equalized_signal = conv(ISI_signal,h);

L = length(equalized_signal);
w = 20*36/(119*5)*(5/6)*((6/5).^[-L:-1]);
equalized_signal = conv(equalized_signal,w);

L = length(equalized_signal);
w = 20/(68*2)*((-1/2).^[0:L-1]);
equalized_signal = conv(equalized_signal,w);

group_delay = 126;
% gain = 1/0.023041;
equalized_signal = equalized_signal(group_delay:end);
equalized_signal = equalized_signal(1:60);
equalized_signal = equalized_signal/max(equalized_signal);
snr_AWGN = SNR(signal,equalized_signal)

figure();
subplot(3,1,1);stem(signal);title('AWGN Transmitted Signal');grid on;
subplot(3,1,2);stem(ISI_signal);title('AWGN Recieced signal');grid on;
subplot(3,1,3);stem(equalized_signal);title('AWGN Recover Causal');grid on;

function isi_sig = ISI_system(signal, multipath, ADD_AWGN, SNR_DB)
  t_sig = trans_branch(signal);
  multipath_effect = multipath;
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
  f_DAC = 5;
  f_DMA = 4;
  srrc_4 =srrc_pulse(4, 5, 1);
  srrc_5 =srrc_pulse(5, 5, 1);
  t_DAC_sig = conv(DAC(sig, f_DAC), srrc_5, 'same');
  t_DMA_sig = conv(DAC(t_DAC_sig, f_DMA), srrc_4, 'same');
  t = [0:length(t_DMA_sig)-1];

  %%%%%%%%%%%%%%
  trans_sig = real(t_DMA_sig .* exp(j*2*pi*fc/fs*t));
  %%%%%%%%%%%%%%
end

function rcv_sig = recieve_branch(demod_sig)
  fc = 16*10^6;
  fs = 64*10^6;
  f_DAC = 5;
  f_DMA = 4;
  srrc_4 =srrc_pulse(4, 5, 1);
  srrc_5 =srrc_pulse(5, 5, 1);

  t = [0:length(demod_sig)-1];
  demod_sig = demod_sig .* exp(-j*2*pi*fc/fs*t);

  f_sig = conv(demod_sig, srrc_4, 'same');
  r_DMA_sig = ADC(f_sig,f_DMA);

  r_DAC_sig = conv(r_DMA_sig,srrc_5, 'same');
  rcv_sig = ADC(r_DAC_sig,f_DAC);
  rcv_sig = real(rcv_sig);
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
  SNR_N = 10^(SNR_DB/10); % SNR enery to linear scale
  SYME = sum(abs(x).^2) / L;
  N0 = SYME / SNR_N;      % Noise spectral Density
  if isreal(x)
    n = sqrt(N0) * randn(1,L);
  else
    n = sqrt(N0/2) * (randn(1,L)+i*randn(1,L));
  end
  y = x + n;
end

function snr = SNR(x,y)
  snr = 10*log10(mean(x.^2) / mean((y-x).^2));
end
