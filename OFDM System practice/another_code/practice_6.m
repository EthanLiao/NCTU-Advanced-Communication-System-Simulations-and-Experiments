clf; clear all; close all;

fb = 1e6;
fd = 1*150e3 /fb;
fIF = 2e6 / fb;
freq_DAC = 16e6;
freq_DMA = 64e6;
f_DAC = freq_DAC / fb;
f_DMA = freq_DMA / freq_DAC;
fb = fb /fb;

N = 60;
sig = randi([0 1],1,N);
sig(sig==0) = -1;

BT = 0.5;
g_filter = Gfilter(BT, f_DAC);
g_delay = (length(g_filter)-1)/2;

load('./filter/CPFSK_IIR')
group_delay = 10;

srrc_4 = rcosine(1,4,'fir/sqrt',0.3,5);
srrc_4_delay = (length(srrc_4)-1)/2;

g_sig = conv(g_filter, DAC(sig, f_DAC));
g_sig = g_sig(g_delay+1:end-g_delay);
g_sig = g_sig / max(g_sig);

for i=1:length(g_sig)
  sum_g_sig(i) =  sum(g_sig(1:i));
end

cp_sig = exp(1j*2*pi*fd/f_DAC*sum_g_sig);
IF_sig = real(cp_sig .* exp(1j*2*pi*fIF/f_DAC*[0:length(cp_sig)-1]));

t_DMA_sig = filter(CPFSK_IIR, DAC(IF_sig, f_DMA));
t_DMA_sig = t_DMA_sig(group_delay:end);
t_DMA_sig = t_DMA_sig / max(t_DMA_sig);

% reciever
r_DMA_sig = filter(CPFSK_IIR, t_DMA_sig);
r_DMA_sig = r_DMA_sig(group_delay:end);
r_DMA_sig = ADC(r_DMA_sig, f_DMA);
r_DMA_sig = r_DMA_sig / max(r_DMA_sig);

r_fIF = r_DMA_sig.*exp(-1j*2*pi*fIF/f_DAC*[0:length(r_DMA_sig)-1]);

r_srrc = conv(r_fIF,srrc_4);
r_srrc = r_srrc(srrc_4_delay+1:end-srrc_4_delay);
r_srrc = r_srrc / max(r_srrc);

r_phase = phase(r_srrc);
r_phase_diff = zeros(1,length(r_phase));
r_phase_diff(1) = r_phase(1);
r_phase_diff(2:end) = r_phase(2:end)-r_phase(1:end-1);

r_DAC_sig = conv(r_phase_diff, g_filter);
r_DAC_sig = r_DAC_sig(g_delay+1:end-g_delay);
r_DAC_sig = ADC(r_DAC_sig, f_DAC);
r_DAC_sig = r_DAC_sig / max(r_DAC_sig);

r_DAC_sig(r_DAC_sig>0) = 0.5;
r_DAC_sig(r_DAC_sig<0) = -0.5;

stem(sig);hold on;stem(r_DAC_sig)
legend("transmitted signal","recieved signal")

figure()
plot(phase(cp_sig));hold on;plot(phase(r_srrc));

figure()
[spec,f] = pwelch(t_DMA_sig,[],[],[],1);
spectrum = 10*log10(spec);
center = max(spectrum);
plot(f,spectrum);
hold on;
spec_vec = [0.5/64:0.001:1/64 1/64:0.001:3/64 3/64:0.001:3.5/64];
mask = [(center-26)*ones(1,8) (center)*ones(1,32) (center-26)*ones(1,8)];
plot(spec_vec, mask);

function u_sig = DAC(sig, u_factor)
  u_sig = zeros(1, length(sig)*u_factor);
  u_sig([1:u_factor:end]) = sig;
end

function d_sig = ADC(sig, d_factor)
  d_sig = zeros(1, length(sig));
  d_sig = sig(1:d_factor:end);
end

function g_filter = Gfilter(BT, M)
  t = [-128:128];
  C = sqrt(2*pi/log(2));
  g_filter =  C * exp(-2*pi^2/log(2)*(BT/M)^2*t.^2);
  g_filter = g_filter / sqrt(sum(g_filter.^2));
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
