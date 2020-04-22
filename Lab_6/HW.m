clf
fc = 6
down = 3
N = 6 % number of sinusoid wave
t = 0:0.01:(N)/(2*pi)
load('A1.mat')
load('Recover_HW_Filter_5')
origin_signal = a1.'
% sound(origin_signal,8000)

[h,w] = freqz(Recover_HW_Filter_5,'whole',2001)
FIR_DC_gain = 10.^(20*log10(abs(h))./20)

down_sampling_signal = ADC(origin_signal,down)
up_sampling_signal = DAC(down_sampling_signal,down)
recover_signal = filter(Recover_HW_Filter_5,up_sampling_signal).*down
sound(recover_signal,8000)

subplot(4,2,1);stem(origin_signal);title('origin signal');grid on;
subplot(4,2,2);stem(abs(fft(origin_signal)));title('origin signal FFT');grid on;
subplot(4,2,3);stem(down_sampling_signal);title('down sampling signal');grid on;
subplot(4,2,4);stem(abs(fft(down_sampling_signal)).*down);title('down sampling signal FFT');grid on;
subplot(4,2,5);stem(up_sampling_signal);title('Up sampling signal');grid on;
subplot(4,2,6);stem(abs(fft(up_sampling_signal)));title('Up sampling signal FFT');grid on;
subplot(4,2,7);stem(recover_signal);title('Recover signal');grid on;
subplot(4,2,8);stem(abs(fft(recover_signal)));title('Recover signal FFT');grid on;
% subplot(3,2,5);stem(down_sampling_signal);title('ADC signal');grid on;

function DAC_sig = DAC(origin_signal,up_factor)
DAC_sig = zeros(1,up_factor*length(origin_signal))
DAC_sig([1:up_factor:length(DAC_sig)]) = origin_signal
end

function ADC_sig = ADC(origin_signal,down_factor)
ADC_sig = origin_signal([1:down_factor:length(origin_signal)])
end
