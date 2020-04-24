function dft_data = dft(sig,dft_size)
  dft_data = fft(sig,dft_size)
  dft_data = ifftshift(dft_data)
end
