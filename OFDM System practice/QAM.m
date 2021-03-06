% This function performs mapping of every 4 data bits into one QAM symbol
function data_modulated = QAM_MAPPER(data)
  for i=1:length(data)
    if data(i) == 0
      data_modulated(i) = -3+3*j
      continue
    elseif data(i) == 1
      data_modulated(i) = -3+1*j
      continue
    elseif data(i) == 2
      data_modulated(i) = -3-1*j
      continue
    elseif data(i) == 3
      data_modulated(i) = -3-3*j
      continue
    elseif data(i) == 4
      data_modulated(i) = -1+3*j
      continue
    elseif data(i) == 5
      data_modulated(i) = -1+1*j
      continue
    elseif data(i) == 6
      data_modulated(i) = -1-3*j
      continue
    elseif data(i) == 7
      data_modulated(i) = -1-1*j
      continue
    elseif data(i) == 8
      data_modulated(i) = 3+3*j
      continue
    elseif data(i) == 9
      data_modulated(i) = 3+1*j
      continue
    elseif data(i) == 10
      data_modulated(i) = 3-3*j
      continue
    elseif data(i) == 11
      data_modulated(i) = 3-1*j
      continue
    elseif data(i) == 12
      data_modulated(i) = 1+3*j
      continue
    elseif data(i) == 13
      data_modulated(i) = 1+1*j
      continue
    elseif data(i) == 14
      data_modulated(i) = 1-3*j
      continue
    else data(i) == 15
      data_modulated(i) = 1-1*j
      continue
    end
  end
end
