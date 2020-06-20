clear all
close all
clc

% % BPSK
N = 256;
x = randi([0 1],1,N);
x((x == 0)) = -1;

L1 = 16;
L2 = 4;
SNR_dB = 20;

%% transmit
% upsampling
L1 = 16;
x_upsampling = zeros(1,N*L1);
x_upsampling(1:L1:end) = x;
% gaussian filter
L1 = 16;
BT = 0.5;
x_gaussian = conv(x_upsampling,gaussfilter(BT,L1),'same');
% sum
fd = 23*150*10^3;
fb = 10^(6)*L1;
Tb = 1/fb;

for n=1:length(x_gaussian)
    x_sum_tmp(n) = sum(x_gaussian(1:n));
end
x_sum = exp(1j*2*pi*fd*Tb*x_sum_tmp);
% upconversion (IF) + real
f_if = 2*10^6;
fs = 16*10^6;
w_if = 2*pi*(f_if/fs);
% t = 1:length(x_sum);
t = 0:length(x_sum)-1;
x_if = real( x_sum.*exp(1j*w_if*t) );
% DAC
L2 = 4;
x_dac = zeros(1,length(x_if)*L2);
x_dac(1:L2:end) = x_if;
% DMA
L2 = 4;
rf = 1;
x_dma = conv(x_dac,srrc(rf,L2,5),'same');

% noise
%SNR_dB = 10;
signal_power = mean(abs(x).^2);
noise_power = signal_power/10^(SNR_dB/10);
noise = sqrt(noise_power)*randn(1,length(x_dma));
y = x_dma;% + noise;
%% receive
% DMA
L2 = 4;
y_dma = conv(y,srrc(rf,L2,5),'same');
% ADC
L2 = 4;
y_adc = y_dma(1:L2:end);
% downconversion (IF)
f_if = 2*10^6;
fs = 16*10^6;
w_if = 2*pi*(f_if/fs);
% t = 1:length(y_adc);
t = 0:length(y_adc)-1;
y_if = y_adc.*exp(-1j*w_if*t);
% srrc filter
L1 = 16;
y_srrc = conv(y_if,srrc(rf,L1,5),'same');
% phase taking
y_phase = unwrap(angle(y_srrc));
% diff
y_diff = zeros(1,length(y_phase));
y_diff(2:end) = y_phase(2:end) - y_phase(1:end-1);
y_diff(1) = y_phase(1);
% gaussian filter
BT = 0.5;
L1 = 16;
y_gaussian = conv(y_diff,gaussfilter(BT,L1),'same');
% downsampling
L1 = 16;
y_downsampling = y_gaussian(1:L1:end);
receive_signal = y_downsampling;

% detect
receive_signal((receive_signal > 0) ) = 1;
receive_signal((receive_signal < 0) ) = -1;

figure();
[px,f] = pwelch(y,[],[],[],1);
plot(f,10*log10(px));

% spectrum mask
hold on;
f_mask = ([-3 -2.5 -2.5:.01:-1.5 -1.5:.01:-1 -1:.01:1 1:.01:1.5 1.5:.01:2.5 2.5 3]+2)/64;
max_mask = max(10*log10(px));
mask = [(max_mask-66) (max_mask-66) (max_mask-46)*ones(1,101) (max_mask-26)*ones(1,51) max_mask*ones(1,201) (max_mask-26)*ones(1,51) (max_mask-46)*ones(1,101) (max_mask-66) (max_mask-66)];
plot(f_mask,mask);
axis([0, 0.1, -inf, inf]);


figure()
plot(phase(x_sum));
hold on
plot(phase(y_srrc));
title("phase comparison");

function y = gaussfilter(BT,M)
    n = -128:128;
    C = sqrt(2*pi/log(2));
    y = C*exp(- (2*(pi)^2/log(2))*(BT/M)^2.*n.^2 );
    y = y/sqrt(sum(abs(y).^2));
end
function y = srrc(rf,M,sp)
    n = [-M*sp:M*sp]+0.00001;
    tmp1 = 4*rf/pi;
    tmp2 = cos((1+rf)*pi*n/M);
    tmp3 = M*sin((1-rf)*pi*n/M)/4/rf./n;
    tmp4 = 1-(4*rf*n/M).^2;
    h_srrc = tmp1*(tmp2+tmp3)./tmp4;
    y = h_srrc/sqrt(sum(abs(h_srrc).^2));
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
