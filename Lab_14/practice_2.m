clf;clear all;close all;
h = [0.1, 0.2, 0.3, 0.4];
h_delay = (length(h)-1)/2;

N = 20;
variance = 1;
noise_sig = sqrt(variance)*randn(1, N);
var_sig = mean((noise_sig-mean(noise_sig)).^2)

filt_sig = conv(noise_sig,h);
filt_sig = filt_sig(h_delay+1:end-h_delay);

histogram(filt_sig);title("Dinamic range of input");

function y = add_awgn_noise(x,SNR_DB)
  L = length(x);
  % calculate symbol energy
  SNR = 10^(SNR_DB/10); % SNR enery to linear scale
  SYME = sum(abs(x).^2) / L;
  N0 = SYME / SNR;      % Noise spectral Density
  if isreal(x)
    n = sqrt(N0) * randn(1,L);
  else
    n = sqrt(N0/2) * (randn(1,L)+i*randn(1,L));
  end
  y = x + n;
end
