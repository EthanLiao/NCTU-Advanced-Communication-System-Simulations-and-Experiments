function desub_sig = subcarrier_demapping(sig,dft_size,used_carrier)
  desub_sig = zeros(used_carrier,1);
  sub_carrier_index=[-used_carrier/2:-1 1:used_carrier/2];
  % a=1:length(data_multiplexed)
  for a=1:used_carrier
    desub_sig(a) = sig(sub_carrier_index(a)+dft_size/2+1)
  end
end
