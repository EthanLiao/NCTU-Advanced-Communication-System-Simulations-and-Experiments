clf;clear all;close all;

h = [0.1, 0.2, 0.3, 0.4];
val = [-1:0.01:1];

% floating point operation
filt_sig = float_operation(val, h);

% fixed point operation
NOB = 3;
q_val = real_ADC(val, NOB, [min(val),max(val)]);
q_h = real_ADC(h, NOB, [min(h),max(h)]);
q_filt_sig = float_operation(q_val, q_h);

q_sqnr = SQNR(filt_sig, q_filt_sig)


function y = float_operation(x,h)
  % tap = zeros(1:length(x)-1);
  z_n = x;
  for i = 1:length(h)
    if i == 1
      y = h(i) * z_n;
    else
      z_n = [0 z_n];
      y = [y 0] + h(i) * z_n;
    end
  end
  % tap = zeros(1:length(x)-1)
  % z_1 = [tap x];
  % y = h(1) * [x tap] + h(2) * z_1;
  %
  % z_2 = [tap, z_1];
  % y =  [y tap] + h(3) * z_2;
  %
  % z_3 = [tap z_2];
  % y = [y tap] + h(4) * z_3;

end

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
