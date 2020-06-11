clf;clear all;close all;

h = [0.1, 0.2, 0.3, 0.4];

N = 1000;
variance = 1;
val = sqrt(variance)*randn(1, N);




% fixed point operation
NOB = 5;
q_val = real_ADC(val, NOB, [min(val),max(val)]);
q_h = real_ADC(h, 2, [min(h),max(h)]);
q_direct_filt_sig = direct_fix_operation(q_val, q_h);
q_tran_filt_sig = tran_fix_operation(q_val, q_h);
q_hybrid_filt_sig = hybrid_fix_operation(q_val, q_h);

% floating point operation
filt_sig = float_operation(val, h);

figure();title("non qunatize vs quantize");
histogram(filt_sig);
hold on;
histogram(q_direct_filt_sig);
hold on;
histogram(q_tran_filt_sig);
hold on;
histogram(q_hybrid_filt_sig);
legend("non quantize"," direct quantize", " trans quantize", "hybrid quantize")

% figure();
% histogram(q_h);
% title("param quantize");
sqnr_direct = SQNR(filt_sig, q_direct_filt_sig)
sqnr_tran = SQNR(filt_sig, q_tran_filt_sig)
sqnr_hybrid = SQNR(filt_sig, q_hybrid_filt_sig)


srrc = srrc_pulse(2,2,0.3);

function y = float_operation(x,h)
  tap = 0;
  z_n = x;
  for i = 1:length(h)
    if i == 1
      y = h(i) * z_n;
    else
      z_n = [tap z_n];
      y = [y tap] + h(i) * z_n;
    end
  end
end

function Sout = direct_fix_operation(x,h)

  figure();
  subplot(2,4,1);
  histogram(x);
  title("Sin quantize");

  Sin = h(1) * x;
  S1 = real_ADC(Sin, 6, [min(Sin),max(Sin)]);
  subplot(2,4,2);
  histogram(S1);
  title("S1 quantize");

  % first stage
  z1 = h(2) * [0 x];
  S2 = real_ADC(z1, 6, [min(z1),max(z1)]);
  subplot(2,4,3);
  histogram(S2);
  title("S2 quantize");

  S5 = [S1 0] + S2;
  S5 = real_ADC(S5, 5, [min(S5),max(S5)]);
  subplot(2,4,4);
  histogram(S5);
  title("S5 quantize");

  % second stage
  z2 = h(3) * [0 0 x];
  S3 = real_ADC(z2, 6, [min(z2),max(z2)]);
  subplot(2,4,5);
  histogram(S3);
  title("S3 quantize");

  S6 = [S5 0] + S3;
  S6 = real_ADC(S6, 5, [min(S6),max(S6)]);
  subplot(2,4,6);
  histogram(S6);
  title("S6 quantize");

  % third stage
  z3 = h(4) * [0 0 0 x];
  S4 = real_ADC(z3, 6, [min(z3),max(z3)]);
  subplot(2,4,7);
  histogram(S4);
  title("S4 quantize");

  Sout = [S6 0] + S4;
  Sout = real_ADC(Sout, 5, [min(Sout),max(Sout)]);
  subplot(2,4,8);
  histogram(Sout);
  title("Sout Direct");
end

function Sout = tran_fix_operation(x,h)

  figure();
  subplot(2,4,1);
  histogram(x);
  title("Sin quantize");

  S1 = h(1) * x;
  S1 = real_ADC(S1, 6, [min(S1),max(S1)]);
  subplot(2,4,2);
  histogram(S1);
  title("S1 quantize");


  z1 = [0 S1];
  S2 = real_ADC(z1, 5, [min(z1),max(z1)]);
  subplot(2,4,3);
  histogram(S2);
  title("S2 quantize");

  S3 = h(2)*x;
  S3 = real_ADC(S3, 5, [min(S3),max(S3)]);
  subplot(2,4,4);
  histogram(S3);
  title("S3 quantize");

  S4 = [S3 0]+S2;
  S4 = real_ADC(S4, 5, [min(S4),max(S4)]);
  subplot(2,4,5);
  histogram(S4);
  title("S4 quantize");

  S5 = h(3)*x;
  S5 = real_ADC(S5, 5, [min(S5),max(S5)]);
  subplot(2,4,6);
  histogram(S5);
  title("S5 quantize");

  % third stage
  z2 = [0 S4];
  S6 = z2+[S5 0 0];
  subplot(2,4,7);
  histogram(S6);
  title("S6 quantize");

  S7 = h(4)*x;
  S7 = real_ADC(S7, 5, [min(S7),max(S7)]);
  subplot(2,4,6);
  histogram(S7);
  title("S7 quantize");

  z3 = [0 S6];
  Sout = z3+[S7 0 0 0];
  subplot(2,4,7);
  histogram(Sout);
  title("Sout transpose");
end


function Sout = hybrid_fix_operation(x,h)

  figure();
  subplot(2,4,1);
  histogram(x);
  title("Sin quantize");

  S1 = h(1) * x;
  S1 = real_ADC(S1, 6, [min(S1),max(S1)]);
  subplot(2,4,2);
  histogram(S1);
  title("S1 quantize");


  S2 = h(2)*[0 x];
  S2 = real_ADC(S2, 5, [min(S2),max(S2)]);
  subplot(2,4,3);
  histogram(S2);
  title("S2 quantize");

  S3 = [S1 0] + S2;
  S3 = real_ADC(S3, 5, [min(S3),max(S3)]);
  subplot(2,4,4);
  histogram(S3);
  title("S3 quantize");

  S4 = h(3)*[0 x];
  S4 = real_ADC(S4, 5, [min(S4),max(S4)]);
  subplot(2,4,5);
  histogram(S4);
  title("S4 quantize");

  S5 = [S4 0]+[0 S3];
  S5 = real_ADC(S5, 5, [min(S5),max(S5)]);
  subplot(2,4,6);
  histogram(S5);
  title("S5 quantize");

  S6 = h(4)*[0 0 x];
  subplot(2,4,7);
  histogram(S6);
  title("S6 quantize");

  Sout = S5+S6;
  Sout = real_ADC(Sout, 5, [min(Sout),max(Sout)]);
  subplot(2,4,8);
  histogram(Sout);
  title("Sout hybrid");

  % for compensation
  Sout = [Sout 0];
end

function Sout = srrc_direct_fix_operation(x,h,srrc)

  figure();
  subplot(2,4,1);
  histogram(x);
  title("Sin quantize");

  Sin = h(1) * x;
  S1 = real_ADC(Sin, 6, [min(Sin),max(Sin)]);
  subplot(2,4,2);
  histogram(S1);
  title("S1 quantize");

  % first stage
  z1 = h(2) * [0 x];
  S2 = real_ADC(z1, 6, [min(z1),max(z1)]);
  subplot(2,4,3);
  histogram(S2);
  title("S2 quantize");

  S5 = [S1 0] + S2;
  S5 = real_ADC(S5, 5, [min(S5),max(S5)]);
  subplot(2,4,4);
  histogram(S5);
  title("S5 quantize");

  % second stage
  z2 = h(3) * [0 0 x];
  S3 = real_ADC(z2, 6, [min(z2),max(z2)]);
  subplot(2,4,5);
  histogram(S3);
  title("S3 quantize");

  S6 = [S5 0] + S3;
  S6 = real_ADC(S6, 5, [min(S6),max(S6)]);
  subplot(2,4,6);
  histogram(S6);
  title("S6 quantize");

  % third stage
  z3 = h(4) * [0 0 0 x];
  S4 = real_ADC(z3, 6, [min(z3),max(z3)]);
  subplot(2,4,7);
  histogram(S4);
  title("S4 quantize");

  Sout = [S6 0] + S4;
  Sout = real_ADC(Sout, 5, [min(Sout),max(Sout)]);
  subplot(2,4,8);
  histogram(Sout);
  title("Sout Direct");
end

% function Sout = tran_fix_operation(x,h)
%
%   figure();
%   subplot(2,4,1);
%   histogram(x);
%   title("Sin quantize");
%
%   S1 = h(1) * x;
%   S1 = real_ADC(S1, 6, [min(S1),max(S1)]);
%   subplot(2,4,2);
%   histogram(S1);
%   title("S1 quantize");
%
%
%   z1 = [0 S1];
%   S2 = real_ADC(z1, 5, [min(z1),max(z1)]);
%   subplot(2,4,3);
%   histogram(S2);
%   title("S2 quantize");
%
%   S3 = h(2)*x;
%   S3 = real_ADC(S3, 5, [min(S3),max(S3)]);
%   subplot(2,4,4);
%   histogram(S3);
%   title("S3 quantize");
%
%   S4 = [S3 0]+S2;
%   S4 = real_ADC(S4, 5, [min(S4),max(S4)]);
%   subplot(2,4,5);
%   histogram(S4);
%   title("S4 quantize");
%
%   S5 = h(3)*x;
%   S5 = real_ADC(S5, 5, [min(S5),max(S5)]);
%   subplot(2,4,6);
%   histogram(S5);
%   title("S5 quantize");
%
%   % third stage
%   z2 = [0 S4];
%   S6 = z2+[S5 0 0];
%   subplot(2,4,7);
%   histogram(S6);
%   title("S6 quantize");
%
%   S7 = h(4)*x;
%   S7 = real_ADC(S7, 5, [min(S7),max(S7)]);
%   subplot(2,4,6);
%   histogram(S7);
%   title("S7 quantize");
%
%   z3 = [0 S6];
%   Sout = z3+[S7 0 0 0];
%   subplot(2,4,7);
%   histogram(Sout);
%   title("Sout transpose");
% end
%
%
% function Sout = hybrid_fix_operation(x,h)
%
%   figure();
%   subplot(2,4,1);
%   histogram(x);
%   title("Sin quantize");
%
%   S1 = h(1) * x;
%   S1 = real_ADC(S1, 6, [min(S1),max(S1)]);
%   subplot(2,4,2);
%   histogram(S1);
%   title("S1 quantize");
%
%
%   S2 = h(2)*[0 x];
%   S2 = real_ADC(S2, 5, [min(S2),max(S2)]);
%   subplot(2,4,3);
%   histogram(S2);
%   title("S2 quantize");
%
%   S3 = [S1 0] + S2;
%   S3 = real_ADC(S3, 5, [min(S3),max(S3)]);
%   subplot(2,4,4);
%   histogram(S3);
%   title("S3 quantize");
%
%   S4 = h(3)*[0 x];
%   S4 = real_ADC(S4, 5, [min(S4),max(S4)]);
%   subplot(2,4,5);
%   histogram(S4);
%   title("S4 quantize");
%
%   S5 = [S4 0]+[0 S3];
%   S5 = real_ADC(S5, 5, [min(S5),max(S5)]);
%   subplot(2,4,6);
%   histogram(S5);
%   title("S5 quantize");
%
%   S6 = h(4)*[0 0 x];
%   subplot(2,4,7);
%   histogram(S6);
%   title("S6 quantize");
%
%   Sout = S5+S6;
%   Sout = real_ADC(Sout, 5, [min(Sout),max(Sout)]);
%   subplot(2,4,8);
%   histogram(Sout);
%   title("Sout hybrid");
%
%   % for compensation
%   Sout = [Sout 0];
% end
%

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
  val(val>range_max) = range_max;
  val(val<range_min) = range_min;
  % quantization
  adc_sig = round((val-range_min)./step)*step;
  adc_sig = range_min + adc_sig;
  adc_sig(adc_sig>range_max-offset) = range_max-step;

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


function [phi, t] = srrc_pulse(T, A, a)
t = [-A*T:A*T] + 10^(-8); % in order to avoid division by zero problems at t=0.
  if (a>0 && a<=1)
     num = cos((1+a)*pi*t/T) + sin((1-a)*pi*t/T) ./ (4*a*t/T);
     denom = 1-(4*a*t./T).^2;
     phi = 4*a/(pi*sqrt(T)) * num ./ denom;
  elseif (a==0)
     phi = 1/(sqrt(T)) * sin(pi*t/T)./(pi*t/T);
  end
end
