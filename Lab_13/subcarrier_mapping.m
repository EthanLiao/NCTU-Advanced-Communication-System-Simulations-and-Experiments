% This function maps the every 300 symbols to 300 subcarriers
function subcarrier_data=subcarrier_mapping(data_multiplexed,dft_size,used_carrier)
  subcarrier_data = zeros(dft_size,1);
  sub_carrier_index=[-used_carrier/2:-1 1:used_carrier/2];
  % a=1:length(data_multiplexed)
  for a=1:used_carrier
    subcarrier_data(sub_carrier_index(a)+dft_size/2+1)=data_multiplexed(a);
  end
end
