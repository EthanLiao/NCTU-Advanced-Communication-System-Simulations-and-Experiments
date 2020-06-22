clf;clear all;close all;

N = 5;
sig = randi([0,1],1,N);
sig((sig==0)) = -1;

% transmit sig
t_sig = trans_branch(sig);

% 1x2 System
H = randn(2,1);
H = sqrt(1/2)*H+1j*sqrt(1/2)*H;
mimo_rcv = H*t_sig;



% recieve signal
mimo_rcv_1 = recieve_branch(mimo_rcv(1,:));
mimo_rcv_2 = recieve_branch(mimo_rcv(2,:));

% SIMO detector
inv_H = 1./H.';
mimo_rcv = [mimo_rcv_1 ; mimo_rcv_2];
rcv_beam = inv_H * mimo_rcv;

% signal detect
rcv_beam = real(rcv_beam)/length(inv_H);
rcv_beam(rcv_beam>0) = 1;
rcv_beam(rcv_beam<0) = -1;

subplot(2,1,1);stem(sig);title("transmission signal");
subplot(2,1,2);stem(rcv_beam);title("recieved signal");

function c_sig = trans_branch(sig)
  % modulatoin part
  fc = 16*10^6;
  fs = 64*10^6;
  f_DAC = 16;
  f_DMA = 4;
  srrc_4 =srrc_pulse(4, 5, 0.3);
  srrc_16 =srrc_pulse(16, 5, 0.3);
  t_DAC_sig = conv(DAC(sig, f_DAC), srrc_16, 'same');
  t_DAC_sig = t_DAC_sig / max(t_DAC_sig);

  t_DMA_sig = conv(DAC(t_DAC_sig, f_DMA), srrc_4, 'same');
  t_DMA_sig = t_DMA_sig / max(t_DMA_sig);

  t = [0:length(t_DMA_sig)-1];
  c_sig = real(t_DMA_sig .* exp(1j*2*pi*fc/fs*t));
end

function r_DAC_sig = recieve_branch(demod_sig)
  fc = 16*10^6;
  fs = 64*10^6;
  f_DAC = 16;
  f_DMA = 4;
  srrc_4 =srrc_pulse(4, 5, 0.3);
  srrc_16 =srrc_pulse(16, 5, 0.3);

  t = [0:length(demod_sig)-1];
  demod_sig = demod_sig .* exp(-1j*2*pi*fc/fs*t);

  f_sig = conv(demod_sig, srrc_4, 'same');
  r_DMA_sig = ADC(f_sig,f_DMA);
  r_DMA_sig = r_DMA_sig / max(r_DMA_sig);

  r_DAC_sig = conv(r_DMA_sig, srrc_16, 'same');
  r_DAC_sig = ADC(r_DAC_sig,f_DAC);
  r_DAC_sig = r_DAC_sig / max(r_DAC_sig);
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
