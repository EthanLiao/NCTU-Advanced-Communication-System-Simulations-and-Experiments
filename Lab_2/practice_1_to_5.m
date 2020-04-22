% practice 1

clf;clear all;
fc = 3              % nature frequency of sinusoidal wave
N = 3               % number of sinusoidal wave
t = 0:0.01:N/fc

signal_freq3 = cos(2*pi*fc*t)
fft_signal_freq3 = abs(fft(signal_freq3 ))

fc = 5
N = 3
t = 0:0.01:N/fc

signal_freq5 = cos(2*pi*fc*t)
fft_signal_freq5 = abs(fft(signal_freq5 ))

% practice 2 : use the frequency = 5 signal
N_factor = 2^2
signal_factor = cos(2*pi*fc*t/N_factor)
fft_signal_factor = abs(fft(signal_freq5 ))


% practice 3 : padding zeros and observe spectrum
padding_num = 2
padding_signal = pad(signal_freq5,2)
fft_padding_signal = abs(fft(padding_signal))


% practice 4
rec_len = 10
mid = ceil(rec_len/2)
half = 3
rect = zeros(1,rec_len)
rect(mid-half:mid+half) = 1
fft_rect = abs(fft(rect))

% practice 5
tri_wave = conv(rect,rect)
tri_wave_with_noise = add_awgn_noise(tri_wave,15)
fft_tri_wave = abs(fft(tri_wave_with_noise))

% practice 1 plot
subplot(6,2,1);stem(signal_freq3);title('Practice1 : sinusoidal wave frequency 3')
subplot(6,2,2);stem(fft_signal_freq3);title('Practice1 : fft sinusoidal wave')
subplot(6,2,3);stem(signal_freq5);title('Practice1 : sinusoidal wave frequency 5')
subplot(6,2,4);stem(fft_signal_freq5);title('Practice1 : fft sinusoidal wave')

% practice 2 plot
subplot(6,2,5);stem(signal_factor);title('Practice2 : sinusoidal wave frequency 5 with factor')
subplot(6,2,6);stem(fft_signal_factor);title('Practice2 : fft sinusoidal with factor wave')

% practice 3 plot
subplot(6,2,7);stem(padding_signal);title('Practice3 : padding sinusoidal wave frequency 5')
subplot(6,2,8);stem(fft_padding_signal);title('Practice3 : fft padding sinusoidal wave frequency 5')

% practice 4 plot
subplot(6,2,9);stem(rect);title('Practice4 : rectangular wave')
subplot(6,2,10);plot(fft_rect);title('Practice4 : fft rectangular wave')

% practice 5 plot
subplot(6,2,11);plot(tri_wave_with_noise);title('Practice5 : triangular wave')
subplot(6,2,12);plot(fft_tri_wave);title('Practice5 : fft triangular wave')


function pad_array = pad(signal,zero_num)
  pad_array = zeros(1,(length(signal)+2*zero_num))
  pad_array((1+zero_num):(end-zero_num)) = signal(1:end)
end

function y = add_awgn_noise(x,SNR_DB)
  L = length(x)
  % calculate symbol energy
  SNR = 10^(SNR_DB/10) % SNR enery to linear scale
  SYME = sum(abs(x).^2) / L
  N0 = SYME / SNR      % Noise spectral Density
  if isreal(x)
    n = sqrt(N0) * randn(1,L)
  else
    n = (N0/2) * (randn(1,L)+i*randn(1,L))
  end
  y = x + n
end
