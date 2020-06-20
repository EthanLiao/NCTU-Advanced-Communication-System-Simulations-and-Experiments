clf;clear all;close all;

N = 10;
sig_1 = randi([0,1],1,N);
sig_1((sig_1==0)) = -1;

sig_2 = randi([0,1],1,N);
sig_2((sig_2==0)) = -1;

sig = [sig_1 ; sig_2];


% 2x2 System
H = randn(2,2);
H = sqrt(1/2) * H +1j * sqrt(1/2) * H;
mimo_rcv = H*sig;


SNR_DB = 0;
% recieve signal
mimo_rcv_1 = add_awgn_noise(mimo_rcv(1,:), SNR_DB);
mimo_rcv_2 = add_awgn_noise(mimo_rcv(2,:), SNR_DB);


mimo_rcv = [mimo_rcv_1 ; mimo_rcv_2];

% ------ Zero Forcing -------------
inv_H = inv(H);
rcv_beam = inv_H * mimo_rcv;

% signal detect
rcv_beam = real(rcv_beam);
% rcv_beam = rcv_beam / max(rcv_beam,[],'all');
% rcv_beam(rcv_beam>0.5) = 1;
% rcv_beam(rcv_beam<-0.5) = -1;

% symbol error rate
% zf_symbol_err = 1-mean(sig==rcv_beam, 'all')
zf_snr = SNR(sig,rcv_beam)
% Zero forcing
subplot(2,2,1);stem(sig_1);title("ZF transmission signal branch 1");
subplot(2,2,2);stem(sig_2);title("ZF transmission signal branch 2");
subplot(2,2,3);stem(rcv_beam(1,:));title("ZF recieved signal branch 1");
subplot(2,2,4);stem(rcv_beam(2,:));title("ZF recieved signal branch 2");

% ------------------MMSE Process-----------------

% MMSE Dectector
% calculate Signal to noise power ratio : sigm

sigm = 10^(SNR_DB/10);

W = inv(H'*H + 1/sigm * eye(2))*H';
% W = H'*inv(H*H' + 1/sigm * eye(2));
mmse_rcv_beam = W * mimo_rcv;

% signal detect
mmse_rcv_beam = real(mmse_rcv_beam);
% rcv_beam = rcv_beam / max(rcv_beam,[],'all');
% mmse_rcv_beam(rcv_beam>0.3) = 1;
% mmse_rcv_beam(rcv_beam<-0.3) = -1;

% symbol error rate
% mmse_symbol_err = 1-mean(sig==rcv_beam, 'all')
mmse_snr = SNR(sig,mmse_rcv_beam)

% MMSE
figure();
subplot(2,2,1);stem(sig_1);title("MMSE transmission signal branch 1");
subplot(2,2,2);stem(sig_2);title("MMSE transmission signal branch 2");
subplot(2,2,3);stem(mmse_rcv_beam(1,:));title("MMSE recieved signal branch 1");
subplot(2,2,4);stem(mmse_rcv_beam(2,:));title("MMSE recieved signal branch 2");


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
  SNR_DB = 0;
  t = [0:length(demod_sig)-1];

  % demod_sig = add_awgn_noise(demod_sig, SNR_DB);
  demod_sig = awgn(demod_sig, SNR_DB);
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
  SNR_N = 10^(SNR_DB/10); % SNR enery to linear scale
  % SYME = sum(abs(x).^2) / L;
  SYME = mean(abs(x).^2);
  N0 = SYME / SNR_N;      % Noise spectral Density
  if isreal(x)
    n = sqrt(N0) * randn(1,L);
  else
    n = sqrt(N0/2) * (randn(1,L)+i*randn(1,L));
  end
  y = x + n;
end

function snr = SNR(x,y)
  snr = 10*log10(mean(x.^2, 'all') / mean((y-x).^2, 'all'));
  % snr = 10*log10(var(x, 0 ,'all') / var((y-x), 0 ,'all'));
  % snr = var(x, 0 ,'all') / var((y-x), 0 ,'all');
end
