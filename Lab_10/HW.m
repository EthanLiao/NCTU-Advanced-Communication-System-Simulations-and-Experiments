clf;clear all; close all;

% Generate qpsk signal
N = 1000;
sig = sign(randn(1,N))+j*sign(randn(1,N));


load('IIR_filter');

fc = 8*10^6;
fs = 32*10^6;
fi = 4*10^6;
f_DAC = 16;
f_DMA = 4;
freq_DAC = 16*10^6;
freq_DMA = 64*10^6;
g = 1.5;
phi = pi /3;
group_delay = 8;
ISREAL = true; % indicator for real signal

srrc_16 = srrc_pulse(16,5,0.3);
srrc_4 = srrc_pulse(16,5,0.3);

% trans_branch(sig,f_DAC,srrc,f_DMA,IIR_filter,group_delay)
mod_sig_real = trans_branch(real(sig),f_DAC,srrc_16,f_DMA,IIR_filter,group_delay);
mod_sig_imag = trans_branch(imag(sig),f_DAC,srrc_16,f_DMA,IIR_filter,group_delay);

% modulatoin with carrier;
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


% IF_reciever(sig_with_carrier,fc,fi,f_DMA,f_DAC,freq_DMA,freq_DAC,IIR_filter,srrc_DMA,srrc_DAC)
rcv_sig_real = IF_reciever(sig_with_carrier,ISREAL,group_delay,fc,fi,f_DMA,f_DAC,freq_DMA,freq_DAC,IIR_filter,srrc_16,srrc_16);
rcv_sig_imag = IF_reciever(sig_with_carrier,~ISREAL,group_delay,fc,fi,f_DMA,f_DAC,freq_DMA,freq_DAC,IIR_filter,srrc_16,srrc_16);


% generate signal for comparison
alpha =  1/2*(1+g*exp(j*phi));
beta = 1/2*(1-g*exp(j*phi));
compar_sig = alpha * sig + beta * conj(sig);


rcv_sig = rcv_sig_real+ j*rcv_sig_imag;
original_sig = scatterplot(sig,1,0,'b.');
hold on;
scatterplot(rcv_sig*2.5,1,0,'k*',original_sig);

function trans_sig = trans_branch(sig,f_DAC,srrc,f_DMA,IIR_filter,group_delay)
  % modulatoin part
  t_DAC_sig = conv(DAC(sig,f_DAC),srrc);
  srrc_len = length(srrc);
  t_DAC_sig = t_DAC_sig((srrc_len-1)/2+1:end-(srrc_len-1)/2);
  t_DAC_sig = t_DAC_sig / max(t_DAC_sig);

  trans_sig = filter(IIR_filter,DAC(t_DAC_sig,f_DMA));
  trans_sig = trans_sig(group_delay:end);
  % trans_sig = trans_sig /(max(trans_sig))
end

function rcv_sig = recieve_branch(demod_sig,f_DMA,IIR_filter,f_DAC,srrc,group_delay)
  f_sig = filter(IIR_filter,demod_sig);
  f_sig = f_sig(group_delay:end);
  r_DMA_sig = ADC(f_sig,f_DMA);

  r_DAC_sig = conv(r_DMA_sig,srrc);
  L = length(srrc);
  r_DAC_sig = r_DAC_sig((L-1)/2+1:end-(L-1)/2);
  r_DAC_sig = r_DAC_sig / max(abs(r_DAC_sig));
  rcv_sig = ADC(r_DAC_sig,f_DAC);
end

function rcv_sig = IF_reciever(sig_with_carrier,ISREAL,group_delay,fc,fi,f_DMA,f_DAC,freq_DMA,freq_DAC,IIR_filter,srrc_DMA,srrc_DAC)
  % IF Band
  srrc_DAC_length = length(srrc_DAC);

  t = [0:length(sig_with_carrier)-1];
  carrier = cos(2*pi*(fc-fi)/freq_DMA*t);

  IF_sig = sig_with_carrier .* carrier;
  r_DMA_sig = filter(IIR_filter,IF_sig);
  r_DMA_sig = r_DMA_sig(group_delay:end);
  r_DMA_sig = ADC(r_DMA_sig,f_DMA);

  % demodulatoin with carrier
  t = [0:length(r_DMA_sig)-1];
  phase_shift = exp(1j*2*pi*fi/freq_DAC*t);
  carrier = exp(-1j*2*pi*fi/freq_DAC*t);
  demod_sig = 0;
  if ISREAL
    demod_sig = real(r_DMA_sig .* carrier);
    % demod_sig = r_DMA_sig * sqrt(2) .* cos(2*pi*fi/freq_DAC*t)
  else
    demod_sig = imag(r_DMA_sig .* carrier);
    % demod_sig = r_DMA_sig .* (-sqrt(2) .* sin(2*pi*fi/freq_DAC*t))
  end
  % demodulation
  r_ADC_f_sig = conv(demod_sig,srrc_DAC);
  r_ADC_f_sig = r_ADC_f_sig((srrc_DAC_length-1)/2+1:end-(srrc_DAC_length-1)/2);
  rcv_sig = ADC(r_ADC_f_sig,f_DAC);
  rcv_sig = rcv_sig/max(rcv_sig);

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
