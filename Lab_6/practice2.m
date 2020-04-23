clf
fc = 6
down = 3
up = 6
N = 30 % number of sinusoid wave
t = 0:0.01:(N)/(2*pi)
origin_signal = cos(2*pi*fc*t)
load('./filter/FIR_LP')
load('./filter/IIR_LP')

% up and down sample signal
up_sampling_signal = DAC(origin_signal,up)
down_then_up_signal = DAC(ADC(origin_signal,up),up)

% filter the signal
[h,w] = freqz(FIR_LP,'whole',2001)
FIR_DC_gain = 10.^(20*log10(abs(h))./20)

[h_IIR,w_IIR] = freqz(IIR_LP,'whole',2001)
IIR_DC_gain = 10.^(20*log10(abs(h_IIR))./20)

% Delay the original signal due to phase change by filter
delay_signal = delay(origin_signal,10)

FIR_filter_sig = filter(FIR_LP,up_sampling_signal)
FIR_filter_sig = FIR_filter_sig./max(FIR_filter_sig)
IIR_filter_sig = filter(IIR_LP,up_sampling_signal)
IIR_filter_sig = IIR_filter_sig./max(IIR_filter_sig)

subplot(4,1,1);stem(origin_signal);title('origin signal');grid on;
subplot(4,1,2);stem(FIR_filter_sig);title('FIR filter up sampling signal');grid on;
subplot(4,1,3);stem(IIR_filter_sig);title('IIR filter up sampling signal');grid on;
subplot(4,1,4);stem(delay_signal);title('delay signal');grid on;


function DAC_sig = DAC(origin_signal,up_factor)
DAC_sig = zeros(1,up_factor*length(origin_signal))
DAC_sig([1:up_factor:length(DAC_sig)]) = origin_signal
end

function ADC_sig = ADC(origin_signal,down_factor)
ADC_sig = origin_signal([1:down_factor:length(origin_signal)])
end

function delay_sig = delay(sig,damt)
  delay_sig = [zeros(1,damt) sig]
end
