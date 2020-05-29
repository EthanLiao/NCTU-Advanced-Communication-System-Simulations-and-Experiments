clf;clear all;close all;
val = [-1:0.1:1];
q_val = real_ADC(val, 3, [-1,1]);
q_val_2 = real_ADC(val, 4, [-1,1]);

q_sqnr = SQNR(val, q_val)
q_sqnr_2 = SQNR(val, q_val_2)

SNR_DB = 2;
awgn_q_val = add_awgn_noise(q_val,SNR_DB);

awgn_sqnr = SQNR(val, awgn_q_val)

function adc_sig = real_ADC(val, bits_amt, range)
% val : analog input vector
% bits_amt : total number of bits to quantize input value
% range : quantization input range, i.e. [min_in, max_in]

  q_lev = 2^bits_amt;
  midPnt = q_lev/2;   % center point
  range_max = range(2);
  range_min = range(1);
  step = (range_max-range_min) / q_lev;
  offset = 0.5*step;  % quantize step offset
  % min and max value clamping
  val(find(val>range_max)) = range_max;
  val(find(val<range_min)) = range_min;

  % quantization
  adc_sig = (round((val-range_min)./step))*step;
  adc_sig = range_min + adc_sig;
  adc_sig(find(adc_sig>range_max-offset)) = range_max-step;

end

function sqnr = SQNR(x,x_q)
  sqnr = 10*log10(mean(x.^2)/mean((x-x_q).^2));
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
