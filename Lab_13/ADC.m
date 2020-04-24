function down_data = ADC(sig,DOWN_FACTOR)
  down_data = zeros(1,length(sig))
  down_data = sig(1:DOWN_FACTOR:end)
end
