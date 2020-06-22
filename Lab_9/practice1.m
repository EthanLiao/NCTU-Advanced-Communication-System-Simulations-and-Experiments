clf;clear all; close all;

N = 10;
sig = randi([0,1],1,N);
sig((sig==0)) = -1;

load('./filter/IIR_filter')

symbol_rate = 1e6;
carrier_frequency = 8e6/symbol_rate;
freq_DAC = 16e6/symbol_rate;
freq_DMA = 32e6/symbol_rate;
group_delay = 1;

srrc_16 = srrc_pulse(16, 5, 0.3);
srrc_2 = srrc_pulse(2,5,0.3);

srrc_16_delay = (length(srrc_16)-1)/2;
srrc_2_delay = (length(srrc_2)-1)/2;

% modulatoin part
t_DAC_sig = conv(DAC(sig,freq_DAC),srrc_16);
t_DAC_sig = t_DAC_sig(srrc_16_delay+1:end-srrc_16_delay);
t_DAC_sig = t_DAC_sig/max(t_DAC_sig);


t_DMA_sig = conv(DAC(t_DAC_sig,freq_DMA/freq_DAC),srrc_2);
t_DMA_sig = t_DMA_sig(srrc_2_delay:end-srrc_2_delay);
t_DMA_sig = t_DMA_sig/max(t_DMA_sig);


t_sig = filter(IIR_filter,t_DMA_sig);
t_sig = t_sig(group_delay:end);

% modulatoin with carrier
t = [0:length(t_sig)-1];
carrier = sqrt(2)*exp(1j*2*pi*carrier_frequency/freq_DMA*t);
sig_with_carrier = real(t_sig .* carrier);

% demodulation
demod_sig = real(sig_with_carrier .* conj(carrier));
r_ADC_sig = filter(IIR_filter,demod_sig);
r_ADC_sig = r_ADC_sig(group_delay:end);
r_ADC_sig = ADC(r_ADC_sig,freq_DMA);

rcv_sig = r_ADC_sig./ max(abs(r_ADC_sig));

% detection
r_sig = (rcv_sig>0)+(-(rcv_sig<0));


subplot(2,1,1);stem(sig);title('BPSK Signal');grid on;
subplot(2,1,2);stem(r_sig);title('Recieved Signal');grid on;


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
    y = 4*a/pi * num./denom;
  else
    y = 1/T * sin(pi.*t./T) ./ (pi*t./T);
  end
end
