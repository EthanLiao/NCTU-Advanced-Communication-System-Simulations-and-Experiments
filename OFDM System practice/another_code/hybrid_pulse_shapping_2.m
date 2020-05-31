clf;clear all; close all;

% Generate qpsk signal
N = 1000;
sig = sign(randn(1,N))+j*sign(randn(1,N));
% for observation
sig_real = real(sig);
sig_imag = imag(sig);

load('./filter/IIR_filter');

fc = 0.5*10^6;
fs = 1*10^6;
f_DAC = 16;
f_DMA = 8;
group_delay = 23;

srrc_16 = srrc_pulse(16,5,0.3);

mod_sig_real = trans_branch(real(sig),f_DAC,srrc_16,f_DMA,IIR_filter,group_delay);
mod_sig_imag = trans_branch(imag(sig),f_DAC,srrc_16,f_DMA,IIR_filter,group_delay);



% modulatoin with carrier
t = [0:length(mod_sig_real)-1];
carrier_cos = sqrt(2)*cos(2*pi*fc/fs*t);
carrier_sin = sqrt(2)*sin(2*pi*fc/fs*t);

sig_with_carrier_real = mod_sig_real.*carrier_cos;
sig_with_carrier_imag = mod_sig_imag.*(-carrier_sin);

sig_with_carrier = sig_with_carrier_real + sig_with_carrier_imag;

% demodulation
demod_sig_real = sig_with_carrier .* carrier_cos;
demod_sig_imag = sig_with_carrier .* (-carrier_sin);

% recieve signal
rcv_sig_real = recieve_branch(demod_sig_real,f_DMA,IIR_filter,f_DAC,srrc_16,group_delay);
rcv_sig_real = (rcv_sig_real>0)-(rcv_sig_real<0);

rcv_sig_imag = recieve_branch(demod_sig_imag,f_DMA,IIR_filter,f_DAC,srrc_16,group_delay);
rcv_sig_imag = (rcv_sig_imag>0)-(rcv_sig_imag<0);

symr = mean(rcv_sig_real == real(sig))
symi = mean(rcv_sig_imag == real(sig))
% snr = SNR(sig, rcv_sig_real+1j*rcv_sig_imag)

function trans_sig = trans_branch(sig,f_DAC,srrc,f_DMA,IIR_filter,group_delay)
  % modulatoin part
  t_DAC_sig = conv(DAC(sig,f_DAC),srrc);
  srrc_delay = (length(srrc)-1)/2;
  t_DAC_sig = t_DAC_sig(srrc_delay+1:end-srrc_delay);
  t_DAC_sig = t_DAC_sig / max(t_DAC_sig);

  trans_sig = filter(IIR_filter,DAC(t_DAC_sig,f_DMA));
  trans_sig = trans_sig(group_delay:end);
end

function rcv_sig = recieve_branch(demod_sig,f_DMA,IIR_filter,f_DAC,srrc,group_delay)
  f_sig = filter(IIR_filter,demod_sig);
  f_sig = f_sig(group_delay:end);
  r_DMA_sig = ADC(f_sig,f_DMA);

  r_DAC_sig = conv(r_DMA_sig,srrc);
  srrc_delay = (length(srrc)-1)/2;
  r_DAC_sig = r_DAC_sig(srrc_delay+1:end-srrc_delay);
  r_DAC_sig = r_DAC_sig / max(abs(r_DAC_sig));
  rcv_sig = ADC(r_DAC_sig,f_DAC);
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

function snr = SNR(tr_sig,rcv_sig)
  snr = mean(rcv_sig.^2)/mean((rcv_sig-tr_sig).^2)
end
