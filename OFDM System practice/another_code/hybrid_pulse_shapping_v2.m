clf;clear all; close all;
% Generate qpsk signal
N = 20;
sig = randi([0,1],1,N);
sig((sig==0)) = -1;


fc = 0.25e6;
fs = 1e6;
f_DAC = 16;
f_DMA = 8;

load('./filter/IIR_filter');
group_delay = 10; % gorup delay of the IIR filter

srrc_16 = srrc_pulse(16,5,0.3);
srrc_8 = srrc_pulse(8,5,0.3);

% transmit branch
mod_sig_real = trans_branch(sig,f_DAC, srrc_16, f_DMA, srrc_8, IIR_filter, group_delay);

% recieve signal
rcv_sig = recieve_branch(mod_sig_real, f_DMA,IIR_filter, srrc_8, f_DAC,srrc_16, group_delay);

snr = SNR(sig, rcv_sig)



subplot(2,1,1);stem(sig);title('BPSK Signal');grid on;
subplot(2,1,2);stem(rcv_sig);title('Recieved Signal');grid on;


function trans_sig = trans_branch(sig,f_DAC,srrc,f_DMA,srrc_DMA,IIR_filter,group_delay)
  % modulatoin part
  srrc_delay = (length(srrc)-1)/2;
  t_DAC_sig = conv(DAC(sig,f_DAC),srrc);
  t_DAC_sig = t_DAC_sig(srrc_delay+1:end-srrc_delay);
  t_DAC_sig = t_DAC_sig/max(t_DAC_sig);

  % srrc_DMA_delay = (length(srrc_DMA)-1)/2;
  % t_DMA_sig = conv(DAC(t_DAC_sig,f_DMA), srrc_DMA);
  % t_DMA_sig = t_DMA_sig(srrc_DMA_delay+1:end-srrc_DMA_delay);
  % t_DMA_sig = t_DMA_sig/max(t_DMA_sig);
  %
  % trans_sig = filter(IIR_filter,t_DMA_sig);
  % trans_sig = trans_sig(group_delay:end);
  % % trans_sig = trans_sig / max(trans_sig);


  trans_sig = filter(IIR_filter,t_DAC_sig);
  trans_sig = trans_sig(group_delay:end);
  % trans_sig = trans_sig / max(trans_sig);
end

function rcv_sig = recieve_branch(demod_sig,f_DMA,IIR_filter,srrc_DMA,f_DAC,srrc,group_delay)
  f_sig = filter(IIR_filter,demod_sig);
  f_sig = f_sig(group_delay:end);
  % f_sig = max(f_sig);

  % srrc_DMA_delay = (length(srrc_DMA)-1)/2;
  % r_DMA_sig = conv(f_sig, srrc_DMA);
  % r_DMA_sig = r_DMA_sig(srrc_DMA_delay+1:end-srrc_DMA_delay);
  % r_DMA_sig = r_DMA_sig /max(r_DMA_sig);
  % r_DMA_sig = ADC(r_DMA_sig,f_DMA);
  %
  % srrc_delay = (length(srrc)-1)/2;
  % r_DAC_sig = conv(r_DMA_sig,srrc);
  % r_DAC_sig = r_DAC_sig(srrc_delay+1:end-srrc_delay);
  % % r_DAC_sig = r_DAC_sig /max(r_DAC_sig);
  % rcv_sig = ADC(r_DAC_sig,f_DAC);
  %
  % % be caution that we should modify gain here
  % rcv_sig = rcv_sig / max(abs(rcv_sig));

  srrc_delay = (length(srrc)-1)/2;
  r_DAC_sig = conv(f_sig,srrc);
  r_DAC_sig = r_DAC_sig(srrc_delay+1:end-srrc_delay);
  % r_DAC_sig = r_DAC_sig /max(r_DAC_sig);
  rcv_sig = ADC(r_DAC_sig,f_DAC);

  % be caution that we should modify gain here
  rcv_sig = rcv_sig / max(abs(rcv_sig));
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
  if isreal(tr_sig)
    snr = 10*log10(mean(rcv_sig.^2)/mean((rcv_sig-tr_sig).^2));
  else
    snr = 10*log10(mean(rcv_sig.^2)/mean(rcv_sig-tr_sig.^2));
  end
end
