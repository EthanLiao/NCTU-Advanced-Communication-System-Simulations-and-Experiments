function y = add_awgn_noise(x,SNR_DB)
  L = length(x)
  snr = 10^(SNR_DB/10)
  SYME = sum(abs(x).^2) / L
  N0 = SYME / snr
  if isreal(x)
    n = sqrt(N0)*randn(1,L)
  else
    n = (N0/2) * (randn(1,L)+j*randn(1,L))
  y = x+n
end
