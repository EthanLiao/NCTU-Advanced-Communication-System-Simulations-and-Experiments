clf;clear all; close all;

% Generate qpsk signal
N = 20;
sig = sign(randn(1,N))+j*sign(randn(1,N));

load('IIR_filter');

fc = 16*10^6;
fs = 64*10^6;
f_DAC = 16;
f_DMA = 4;
g = 6;
phi = 0;
group_delay = 3;

srrc_16 = srrc_pulse(16,5,0.3);
srrc_16_len = length(srrc_16);

mod_sig_real = trans_branch(real(sig),f_DAC,srrc_16,f_DMA,IIR_filter,group_delay);
mod_sig_imag = trans_branch(imag(sig),f_DAC,srrc_16,f_DMA,IIR_filter,group_delay);

% modulatoin with carrier
t = [0:length(mod_sig_real)-1];
carrier_cos = sqrt(2)*cos(2*pi*fc/fs*t);
carrier_sin = sqrt(2)*sin(2*pi*fc/fs*t);
carrier_sin_im = sqrt(2)*sin(2*pi*fc/fs*t+phi)*g;


% compensation
aIE = mod_sig_real;
aQE = mod_sig_imag;

imbalance_signal = [aIE;aQE];
H = [1 -g*sin(phi) ; 0 g*cos(phi)];
compensate_sig = inv(H) * imbalance_signal;

% conduct a demodulation
im_sig_real = compensate_sig(1,:);
im_sig_imag = compensate_sig(2,:);
sig_with_carrier = im_sig_real.*carrier_cos+im_sig_imag.*(-carrier_sin);

% CFO
t = [0:length(sig_with_carrier)-1];
cfo = 0.01*fc;
carrier_cos_cfo = sqrt(2)*cos(2*pi*(fc+cfo)/fs*t);
carrier_sin_cfo = sqrt(2)*sin(2*pi*(fc+cfo)/fs*t);

% demodulation
demod_sig_real = sig_with_carrier .* carrier_cos_cfo;
demod_sig_imag = sig_with_carrier .* (-carrier_sin_cfo);

rcv_sig_real = recieve_branch(demod_sig_real,f_DMA,IIR_filter,f_DAC,srrc_16,group_delay);
rcv_sig_imag = recieve_branch(demod_sig_imag,f_DMA,IIR_filter,f_DAC,srrc_16,group_delay);
scatterplot(rcv_sig_real+1j*rcv_sig_imag);

% correlate the incorrect signal

rcv_sig = rcv_sig_real + 1j*rcv_sig_imag;

% cfo_rcv_sig = rcv_sig;
% cfo_rcv_sig(real(cfo_rcv_sig)>0 & imag(cfo_rcv_sig)>0) = 1+1j;
% cfo_rcv_sig(real(cfo_rcv_sig)>0 & imag(cfo_rcv_sig)<0) = 1-1j;
% cfo_rcv_sig(real(cfo_rcv_sig)<0 & imag(cfo_rcv_sig)>0) = -1+1j;
% cfo_rcv_sig(real(cfo_rcv_sig)<0 & imag(cfo_rcv_sig)<0) = -1-1j;

[cfo_per,cfo_db] = EVM(rcv_sig,sig);
rcv_sig = sqrt(2)*rcv_sig.*exp(1j*2*pi*cfo_per*fc/(10^6)*[0:length(rcv_sig)-1]);

figure()
subplot(2,1,1);stem(real(sig));hold on;stem(rcv_sig_real);title("non comp real signal")
subplot(2,1,2);stem(imag(sig));hold on;stem(rcv_sig_imag);title("non comp imag signal")


figure();
subplot(2,1,1);stem(real(sig));hold on;stem(real(rcv_sig)/max(real(rcv_sig)));title("comp real signal")
subplot(2,1,2);stem(imag(sig));hold on;stem(imag(rcv_sig)/max(imag(rcv_sig)));title("comp imag signal")



function trans_sig = trans_branch(sig,f_DAC,srrc,f_DMA,IIR_filter,group_delay)
  % modulatoin part
  srrc_delay = (length(srrc)-1)/2;
  t_DAC_sig = conv(DAC(sig,f_DAC),srrc);
  t_DAC_sig = t_DAC_sig(srrc_delay+1:end-srrc_delay);
  t_DAC_sig = t_DAC_sig / max(t_DAC_sig);

  trans_sig = filter(IIR_filter,DAC(t_DAC_sig,f_DMA));
  trans_sig = trans_sig(group_delay:end);
  trans_sig = trans_sig /max(trans_sig);
end

function rcv_sig = recieve_branch(demod_sig,f_DMA,IIR_filter,f_DAC,srrc,group_delay)

  f_sig = filter(IIR_filter,demod_sig);
  f_sig = f_sig(group_delay:end);
  r_DMA_sig = ADC(f_sig,f_DMA);
  r_DMA_sig = r_DMA_sig / max(r_DMA_sig);

  srrc_delay = (length(srrc)-1)/2;
  r_DAC_sig = conv(r_DMA_sig,srrc);
  r_DAC_sig = r_DAC_sig(srrc_delay+1:end- srrc_delay);
  rcv_sig = ADC(r_DAC_sig,f_DAC);
  rcv_sig = rcv_sig / max(rcv_sig);
end

function DAC_sig = DAC(origin_sig,up_factor)
  DAC_sig = zeros(1,up_factor*length(origin_sig));
  DAC_sig(1:up_factor:end) = origin_sig;
end

function ADC_sig = ADC(origin_sig,down_factor)
  ADC_sig = zeros(1,length(origin_sig));
  ADC_sig = origin_sig(1:down_factor:end);
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

function [EVM_per,EVM_dB] = EVM(real_sig, ideal_sig)
  EVM_per = sqrt(mean(abs(ideal_sig-real_sig).^2)/mean(abs(ideal_sig).^2))/1000;
  % EVM_per = angle(mean(ideal_sig-real_sig));
  EVM_dB = 10*log10(mean(abs(ideal_sig-real_sig).^2)/mean(abs(ideal_sig).^2));
end
