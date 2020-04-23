fc = 6
down = 10             % down sampling factor
up = 3
N = 6                 % number of sinusoid wave
t = 0:0.01:(N)/(2*pi)

origin_signal = cos(2*pi*fc*t)
down_sampling_signal = ADC(origin_signal,down)
up_sampling_signal = DAC(origin_signal,up)

subplot(3,2,1);stem(origin_signal);title('origin signal');grid on;
subplot(3,2,3);stem(down_sampling_signal);title('down sampling signal');grid on;
subplot(3,2,5);stem(up_sampling_signal);title('up sampling signal');grid on;

subplot(3,2,2);stem(abs(fft(origin_signal)));title('origin signal');grid on;
subplot(3,2,4);stem(abs(fft(down_sampling_signal)));title('down sampling signal');grid on;
subplot(3,2,6);stem(abs(fft(up_sampling_signal)));title('up sampling signal');grid on;

function down_sig = ADC(signal,down)
  down_sig = zeros(1,length(signal))
  down_sig = signal(1:down:end)
end

function up_sig = DAC(signal,up)
  up_sig = zeros(1,up*length(signal))
  up_sig(1:up:end) = signal
end
