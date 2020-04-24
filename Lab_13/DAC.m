function d_sig = DAC(sig,up_factor)
  d_sig = zeros(1,length(sig)*up_factor)
  d_sig(1:up_factor:end) = sig
end
