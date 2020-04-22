clf
fc = 6
down = 32
up = 32
N = 6 % number of sinusoid wave
t = 0:0.01:(N)/(2*pi)
origin_signal = cos(2*pi*fc*t)
down_sampling_signal = origin_signal([1:down:length(origin_signal)])
up_sampling_signal = zeros(1,up*length(origin_signal))
up_sampling_signal([1:up:length(up_sampling_signal)]) = origin_signal

DAC_orig_signal = DAC(origin_signal,up)

[h,w] = freqz(FIR_LP,'whole',2001)
FIR_DC_gain = 10.^(20*log10(abs(h))./20)


subplot(3,2,1);stem(origin_signal);title('origin signal');grid on;
subplot(3,2,2);stem(abs(fft(origin_signal)));title('origin signal');grid on;
subplot(3,2,3);stem(ADC(32.*filter(FIR_LP,DAC_orig_signal),down));title('DMA_filter');grid on;
subplot(3,2,4);stem(abs(fft(ADC(32.*filter(FIR_LP,DAC_orig_signal),down))));title('DMA_filter');grid on;
subplot(3,2,5);stem(down_sampling_signal);title('ADC signal');grid on;

function DAC_sig = DAC(origin_signal,up_factor)
DAC_sig = zeros(1,up_factor*length(origin_signal))
DAC_sig([1:up_factor:length(DAC_sig)]) = origin_signal
end

function ADC_sig = ADC(origin_signal,down_factor)
ADC_sig = origin_signal([1:down_factor:length(origin_signal)])
end
