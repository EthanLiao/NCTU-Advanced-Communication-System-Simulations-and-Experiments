clf;clear all

fc = 6
down = 32
up = 32
N = 6
t = (0:0.01:N)/(fc)
origin_signal = cos(2*pi*fc*t)
load('./filter/FIR_LP')

% up then filter it then down sample it

down_sampling_signal = ADC(origin_signal,down)
DAC_orig_signal = DAC(origin_signal,up)
recover_signal = up.*ADC(filter(FIR_LP,DAC_orig_signal),down)

% [h,w] = freqz(FIR_LP,'whole',2001)
% FIR_DC_gain = 10.^(20*log10(abs(h))./20)


subplot(3,2,1);stem(origin_signal);title('origin signal');grid on;
subplot(3,2,2);stem(abs(fft(origin_signal)));title('origin signal');grid on;
subplot(3,2,3);stem(recover_signal);title('recover signal by DMA filter');grid on;
subplot(3,2,4);stem(abs(fft(recover_signal)));title('DMA_filter');grid on;
subplot(3,2,5);stem(down_sampling_signal);title('ADC signal');grid on;

function DAC_sig = DAC(origin_signal,up_factor)
DAC_sig = zeros(1,up_factor*length(origin_signal))
DAC_sig([1:up_factor:length(DAC_sig)]) = origin_signal
end

function ADC_sig = ADC(origin_signal,down_factor)
ADC_sig = origin_signal([1:down_factor:length(origin_signal)])
end
