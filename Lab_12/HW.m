% Parametors
clf;close all;close all;
% notice that fs fd fIF should be normalized
fs = 1*10^6;
fc = 8*10^6/fs;
fd = 1*150*10^3/fs;
fIF = 2*10^6/fs;
f_DAC = 16*10^6/fs;
f_DMA = 64*10^6/fs;
fs = fs/fs;
Tb = 1/fs;
BT = 0.5;

N = 13;
signal = randi([0 1],1,N);
signal(signal==0) = -1;

gau_filter = Gfilter(BT,f_DAC,Tb);
gau_delay = (length(gau_filter)-1)/2;
figure()
plot(gau_filter);title("Gaussian filter");

srrc_16 = srrc_pulse(f_DAC,10,1);
srrc_16_delay = (length(srrc_16)-1)/2;

srrc_4 = srrc_pulse(f_DMA/f_DAC,10,1);
srrc_4_delay = (length(srrc_4)-1)/2;

% Gaussian filter
t_DAC_sig = DAC(signal, f_DAC/fs);
g_sig = conv(t_DAC_sig, gau_filter);
g_sig = g_sig(gau_delay+1:end-gau_delay);
% g_sig = g_sig / max(g_sig);

% sum
for i = 1:length(g_sig)
  sum_gsig(i) = sum(g_sig(1:i));
end
cp_sig = exp(1j*2*pi*fd*1/fs*sum_gsig);



% IF modulation
t = [0:length(cp_sig)-1];
IF_sig = real(cp_sig.*exp(1j*2*pi*fIF/(f_DAC/fs)*t));


% DAC/F
t_DMA_sig = conv(DAC(IF_sig, f_DMA/f_DAC), srrc_4);
t_DMA_sig = t_DMA_sig(srrc_4_delay+1:end-srrc_4_delay);
t_DMA_sig = t_DMA_sig/max(t_DMA_sig);

% modulation with Carrier
t = [0:length(t_DMA_sig)-1];
carrier = exp(1j*2*pi*(fc-fIF)/fs*t);
carrier_sig = t_DMA_sig .*carrier;

% AWGN
SNR_DB = 1;
noise_sig = AWGN(carrier_sig,SNR_DB);

% demodulation
t = [0:length(noise_sig)-1];
carrier = exp(-1j*2*pi*(fc-fIF)/fs*t);
demod_sig = noise_sig.*carrier;

% F/ADC
r_DMA_sig = conv(t_DMA_sig, srrc_4);
r_DMA_sig = r_DMA_sig(srrc_4_delay+1:end-srrc_4_delay);
r_DMA_sig = r_DMA_sig / max(r_DMA_sig);
r_DMA_sig = ADC(r_DMA_sig, f_DMA/f_DAC);


% IF demodulation
t = [0:length(r_DMA_sig)-1];
dmod_IF_sig = r_DMA_sig .* exp(-1j*2*pi*fIF/(f_DAC/fs)*t);

% F/PH
re_srrc = conv(dmod_IF_sig,srrc_4);
re_srrc = re_srrc(srrc_4_delay+1:end-srrc_4_delay);
re_srrc = re_srrc/max(re_srrc);


% phase and diff
re_phase = phase(re_srrc);
re_phase_diff(1) = re_phase(1)/(2*pi*fd*1/fs);
for i = 2:length(re_phase)
  re_phase_diff(i) = (re_phase(i)-re_phase(i-1))/(2*pi*fd*1/fs);
end

% Gaussian Filter and ADC
r_sig = conv(re_phase_diff, gau_filter);
r_sig = r_sig(gau_delay+1:end-gau_delay);
r_sig = r_sig / max(r_sig);

r_sig = ADC(r_sig, f_DAC/fs);
r_sig(r_sig<0) = -0.5;
r_sig(r_sig>0) = 0.5;

figure()
stem(signal)
hold on
stem(r_sig);title('Recieved signal');
legend('transmitted','recieved');


figure()
plot(phase(cp_sig));
hold on
plot(phase(re_srrc));
title("phase comparison");

figure()
plot(abs(fft(cp_sig)));
hold on
plot(abs(fft(re_srrc)));
title('frequency domain of gaussian modulated and demodulated signal');
legend("modulated", "demodulated");

function y = AWGN(x,N_dB)
  L = length(x);
  SNR_DB = 10^(N_dB/10);
  SYME = sum(abs(x).^2)/L;
  N0 = SYME / SNR_DB;
  if isreal(x)
    n = sqrt(N0) * randn(1,L);
  else
    n = sqrt(N0/2) * (randn(1,L)+1j*randn(1,L));
  end
  y = x + n ;
end


% some function for reuse
function up_sig = DAC(sig,up_factor)
  up_sig = zeros(1,length(sig)*up_factor);
  up_sig(1:up_factor:end) = sig;
end

function down_sig = ADC(sig, down_factor)
  down_sig = zeros(1, length(sig));
  down_sig = sig(1:down_factor:end);
end

% Gaussian Filter
function g_filter = Gfilter(BT, M, Tb)
  t = [-64:64];
  B = BT/Tb;
  C = sqrt(2*pi/log(2)) * B;
  g_filter = C*exp(-2*(pi^2)/log(2)*(BT/M)^2*(t.^2));
  g_filter = g_filter ./ sqrt(sum(g_filter.^2));
end


function [y,t] = srrc_pulse(T,A,a)
  t = [-A*T:A*T] + 10^(-8);
  if (a>0 && a<=1)
    num = cos((1+a)*pi*t/T) + T*sin((1-a)*pi*t/T)./(4*a*t);
    denom = 1-(4*a*t/T).^2;
    y = (4*a/pi) * num./denom;
  else
    y = 1/T * sin(pi.*t./T) ./ (pi*t./T);
  end
end

function PHI=phase(G)
  %PHASE  Computes the phase of a complex vector
  %
  %   PHI=phase(G)
  %
  %   G is a complex-valued row vector and PHI is returned as its
  %   phase (in radians), with an effort made to keep it continuous
  %   over the pi-borders.

  %   L. Ljung 10-2-86
  %   Copyright 1986-2004 The MathWorks, Inc.
  %   $Revision: 1.5.4.2 $  $Date: 2004/07/31 23:24:49 $

  %PHI = unwrap(angle(G));
  [nr,nc] = size(G);
  if min(nr,nc) > 1
      error(sprintf(['PHASE applies only to row or column vectors.'...
          '\nFor matrices you have to decide along which dimension the'...
          '\nphase should be continuous.']))
  end
  if nr>nc
      G = G.';
  end
  PHI=atan2(imag(G),real(G));
  N=length(PHI);
  DF=PHI(1:N-1)-PHI(2:N);
  I=find(abs(DF)>3.5);
  for i=I
      if i~=0,
          PHI=PHI+2*pi*sign(DF(i))*[zeros(1,i) ones(1,N-i)];
      end
  end
  if nr>nc
      PHI = PHI.';
  end
end
