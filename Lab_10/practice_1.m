clf;clear all; close all;

% Generate qpsk signal
N = 10;
sig = sign(randn(1,N))+j*sign(randn(1,N));
load('IIR_filter');

fc = 16*10^6;
fs = 64*10^6;
f_DAC = 16;
f_DMA = 4;
g = 1;
phi = pi/3;
group_delay = 6;

srrc_16 = srrc_pulse(16,5,0.3);
srrc_16_len = length(srrc_16);

mod_sig_real = trans_branch(real(sig), f_DAC, srrc_16, f_DMA,IIR_filter, group_delay);
mod_sig_imag = trans_branch(imag(sig), f_DAC, srrc_16, f_DMA,IIR_filter, group_delay);

% modulatoin with carrier
t = [0:length(mod_sig_real)-1];
carrier_cos = sqrt(2)*cos(2*pi*fc/fs*t);
carrier_sin = sqrt(2)*sin(2*pi*fc/fs*t);
carrier_sin_im = sqrt(2)*sin(2*pi*fc/fs*t+phi)*g;

% plot(carrier_sin_im)
sig_with_carrier_real = mod_sig_real.*carrier_cos;
sig_with_carrier_imag = mod_sig_imag.*(-carrier_sin_im);

sig_with_carrier = sig_with_carrier_real + sig_with_carrier_imag;

% demodulation
demod_sig_real = sig_with_carrier .* carrier_cos;
demod_sig_imag = sig_with_carrier .* (-carrier_sin);

rcv_sig_real = recieve_branch(demod_sig_real,f_DMA,IIR_filter,f_DAC,srrc_16,group_delay);
rcv_sig_imag = recieve_branch(demod_sig_imag,f_DMA,IIR_filter,f_DAC,srrc_16,group_delay);


% generate signal for comparison
alpha =  1/2*(1+g*exp(j*phi));
beta = 1/2*(1-g*exp(j*phi));
compar_sig = alpha * sig + beta * conj(sig);

subplot(2,1,1);stem(real(sig));hold on;stem(rcv_sig_real);hold on;stem(real(compar_sig));title('QPSK Real Signal');grid on;legend('original signal','imbalance signal','generated siganl');
subplot(2,1,2);stem(imag(sig));hold on;stem(rcv_sig_imag);hold on;stem(imag(compar_sig));title('QPSK Image Signal');grid on;legend('original signal','imbalance signal','generated siganl');

rcv_sig = rcv_sig_real+ j*rcv_sig_imag;
original_sig = scatterplot(sig,1,0,'b.');
hold on;
scatterplot(rcv_sig*2.5,1,0,'k*',original_sig)

function trans_sig = trans_branch(sig,f_DAC,srrc,f_DMA,IIR_filter,group_delay)
  % modulatoin part
  srrc_delay = (length(srrc)-1)/2;

  t_DAC_sig = conv(DAC(sig,f_DAC),srrc);
  t_DAC_sig = t_DAC_sig(srrc_delay+1:end-srrc_delay);
  t_DAC_sig = t_DAC_sig / max(t_DAC_sig);

  trans_sig = filter(IIR_filter,DAC(t_DAC_sig,f_DMA));
  trans_sig = trans_sig(group_delay:end);
  trans_sig = trans_sig /(max(trans_sig));
end

function rcv_sig = recieve_branch(demod_sig,f_DMA,IIR_filter,f_DAC,srrc,group_delay)

  f_sig = filter(IIR_filter,demod_sig);
  f_sig = f_sig(group_delay:end);
  f_sig = ADC(f_sig,f_DMA);
  f_sig = f_sig / max(f_sig);

  srrc_delay = (length(srrc)-1)/2;
  r_DAC_sig = conv(f_sig,srrc);
  r_DAC_sig = r_DAC_sig(srrc_delay+1:end-srrc_delay);
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


function rcv_sig = IF_reciever(sig_with_carrier,ISREAL,fc,fi,f_DMA,f_DAC,freq_DMA,freq_DAC,srrc_DMA,srrc_DAC)
  % IF Band
  srrc_DMA_length = length(srrc_DMA);
  srrc_DAC_length = length(srrc_DAC);

  t = [0:length(sig_with_carrier)-1];
  carrier = cos(2*pi*(fc-fi)/freq_DMA*t);

  IF_sig = sig_with_carrier .* carrier;
  r_DMA_f_sig = conv(IF_sig,srrc_DMA);
  r_DMA_f_sig = r_DMA_f_sig((srrc_DMA_length-1)/2+1:end-(srrc_DMA_length-1)/2);
  r_DMA_sig = ADC(r_DMA_f_sig,f_DMA);

  % demodulatoin with carrier
  t = [0:length(r_DMA_sig)-1];
  carrier = exp(-1j*2*pi*fi/freq_DAC*t);
  if ISREAL
    demod_sig = real(r_DMA_sig .* carrier);
  else
    demod_sig = imag(r_DMA_sig .* carrier);
  end
  % demodulation
  r_ADC_f_sig = conv(demod_sig,srrc_DAC);
  r_ADC_f_sig = r_ADC_f_sig((srrc_DAC_length-1)/2+1:end-(srrc_DAC_length-1)/2);
  r_ADC_sig = ADC(r_ADC_f_sig,f_DAC);
  rcv_sig = r_ADC_sig./max(abs(r_ADC_sig));

end
