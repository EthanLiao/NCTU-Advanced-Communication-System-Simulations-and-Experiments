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
