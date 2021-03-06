clf;clear all;close all;
h = [0.1, 0.2, 0.3, 0.4];
h_delay = (length(h)-1)/2;

N = 1000;
variance = 1;
noise_sig = sqrt(variance)*randn(1, N);
noise_sig = [-1:0.01:1];
% noise_sig = randn(1,0, 1);
% var_sig = mean((noise_sig-mean(noise_sig)).^2);

filt_sig = float_operation(noise_sig,h);
histogram(filt_sig);title("Dinamic range of input");


function y = float_operation(x,h)
  % tap = zeros(1:length(x)-1);
  tap = 0;
  z_n = x;
  for i = 1:length(h)
    if i == 1
      y = h(i) * z_n;
    else
      z_n = [tap z_n];
      y = [y tap] + h(i) * z_n;
    end
  end
end

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
