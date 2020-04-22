clf; clear all;
% HW 1
fs_orig = 4 * 10^6
[sin_sig,sin_time,fft_sin] = sinusoid_with_fft(fs_orig)
fs_another = 1.5 * 10^6
[sin_sig,sin_time,fft_sin] = sinusoid_with_fft(fs_another)

figure()
subplot(4,1,1);stem(sin_time,sin_sig,'.-');title('Sinusoid Signal');xlabel('Time(seconds)'); ylabel('Amplitude'); grid on;
subplot(4,1,2);stem(abs(fft_sin),'.-b') ; xlabel('Frequency'); ylabel('Magnitude'); grid on;
subplot(4,1,3);stem(sin_time,sin_sig,'.-');title('Sinusoid Signal');xlabel('Time(seconds)'); ylabel('Amplitude'); grid on;
subplot(4,1,4);stem(abs(fft_sin),'.-b') ; xlabel('Frequency'); ylabel('Magnitude'); grid on;



% HW 2
ele_nums = 200
mid = ceil(ele_nums/2)
half = 30
% rectangular pulse
rect = zeros(1, ele_nums)
rect(mid-half : mid+half) = 1

% Convolve
triangular = conv(rect, rect)
convolve_Length = length(rect)+length(rect)-1

% Add noise to the triangular wave
SNR = 15
Second_moment_noise = moment(triangular,2) / (10^(SNR/10))
noise =  Second_moment_noise*randn(1,convolve_Length) / 3
tri_with_noise = noise + triangular
fft_y =  fft(tri_with_noise)

figure()
subplot(5, 1, 3); plot(abs(fft_y)); grid 'on';

% filtering white noise in the frequency domain
for idx = 1:convolve_Length
  if abs(fft_y(idx)) < 2000
    fft_y(idx) = 0
  end
end

% generate the window
window_num = 400
window_mid = ceil(window_num/2)
window_size = 50

freq_window = zeros(1,window_num)
freq_window(window_mid-window_size : window_mid+window_size) = 1
fft_window = fft(freq_window)

for idx = 1:400
  if abs(fft_window(idx)) < 10000000
    fft_window(idx) = 0
  end
end

real_window = ifft(fft_window)

% calculate second moment
orig_noise_sec = moment(noise,2)
new_noise_sec = moment(ifft(fft_y)-triangular,2)
signal_sec = moment(triangular,2)

% SNR
orig_SNR = 10*log(signal_sec/orig_noise_sec)/log(10)
new_SNR = 10*log(signal_sec/new_noise_sec)/log(10)

mean_square_error = immse(tri_with_noise,triangular)

subplot(5, 1, 1); plot(triangular); grid 'on';
subplot(5, 1, 2); plot(tri_with_noise); grid 'on'; title(['SNR: ',num2str(orig_SNR),' dB']);
subplot(5, 1, 4); plot(abs(fft_y)); grid 'on';
subplot(5, 1, 5); plot(ifft(fft_y)); grid 'on'; title(['SNR: ',num2str(new_SNR),' dB;','MSE: ',num2str(mean_square_error)]);


function [sig,sig_time,fft_sig] = sinusoid_with_fft(fs)
  N = 3
  % freqStep = fs / N
  fc = 10^6 % fc = 1 * freqStep
  sig_time = (0:0.01:N) / fc
  sig = cos(2*pi*fc*sig_time)
  fft_sig = fft(sig)
  fft_sig = fftshift(fft_sig)
end
