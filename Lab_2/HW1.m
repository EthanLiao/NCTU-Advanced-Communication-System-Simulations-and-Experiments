fs_orig = 4 * 10^6
[sin_sig,sin_time,fft_sin] = sinusoid_with_fft(fs_orig)
subplot(4,1,1);stem(sin_time,sin_sig,'.-');title('Sinusoid Signal');xlabel('Time(seconds)'); ylabel('Amplitude'); grid on;
subplot(4,1,2);stem(abs(fft_sin),'.-b') ; xlabel('Frequency'); ylabel('Magnitude'); grid on;

fs_another = 1.5 * 10^6
[sin_sig,sin_time,fft_sin] = sinusoid_with_fft(fs_another)
subplot(4,1,3);stem(sin_time,sin_sig,'.-');title('Sinusoid Signal');xlabel('Time(seconds)'); ylabel('Amplitude'); grid on;
subplot(4,1,4);stem(abs(fft_sin),'.-b') ; xlabel('Frequency'); ylabel('Magnitude'); grid on;

function [sig,sig_time,fft_sig] = sinusoid_with_fft(fs)
  N = 64
  % freqStep = fs / N
  fc = 10^6 % fc = 1 * freqStep
  sig_time = (0:N-1) / fs
  sig = cos(2*pi*fc*sig_time)
  fft_sig = fft(sig)
  fft_sig = fftshift(fft_sig)
end
