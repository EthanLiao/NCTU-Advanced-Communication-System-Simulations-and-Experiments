fc = 6
down = 10
up = 3
N = 6 % number of sinusoid wave
t = 0:0.01:(N)/(2*pi)
origin_signal = cos(2*pi*fc*t)
down_sampling_signal = origin_signal([1:down:length(origin_signal)])
up_sampling_signal = zeros(1,up*length(origin_signal))
up_sampling_signal([1:up:length(up_sampling_signal)]) = origin_signal
subplot(3,2,1);stem(origin_signal);title('origin signal');grid on;
subplot(3,2,3);stem(down_sampling_signal);title('down sampling signal');grid on;
subplot(3,2,5);stem(up_sampling_signal);title('up sampling signal');grid on;

subplot(3,2,2);stem(abs(fft(origin_signal)));title('origin signal');grid on;
subplot(3,2,4);stem(abs(fft(down_sampling_signal)));title('down sampling signal');grid on;
subplot(3,2,6);stem(abs(fft(up_sampling_signal)));title('up sampling signal');grid on;
% subplot(3,1,3);plot(filter(LP,down_sampling_signal));title('filtered signal');grid on;
