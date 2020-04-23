clf;clear all

% fc = 6
% N = 6
% t = (0:0.01:N)/(fc)
% origin_signal = cos(2*pi*fc*t)
load('./filter/FIR_LP')
load('./filter/FIR_LP_2')
down = 32
up = 32
% generate signal
len = 50
mid = ceil(len/2)
half = 10
rect = zeros(1,len)
rect(mid-half:mid+half) = 1
tri = conv(rect,rect)

% up then filter it then down sample it
up_signal = DAC(tri,up)
f_signal = filter(FIR_LP_2,up_signal)
recover_signal = ADC(f_signal,down)./7.9.*21


subplot(2,2,1);stem(tri);title('origin signal');grid on;
subplot(2,2,2);stem(abs(fft(tri)));title('origin signal f domain');grid on;
subplot(2,2,3);stem(recover_signal);title('recover signal by DMA filter');grid on;
subplot(2,2,4);stem(abs(fft(recover_signal)));title('DMA_filter');grid on;

function DAC_sig = DAC(origin_signal,up_factor)
DAC_sig = zeros(1,up_factor*length(origin_signal))
DAC_sig([1:up_factor:length(DAC_sig)]) = origin_signal
end

function ADC_sig = ADC(origin_signal,down_factor)
ADC_sig = origin_signal([1:down_factor:length(origin_signal)])
end
