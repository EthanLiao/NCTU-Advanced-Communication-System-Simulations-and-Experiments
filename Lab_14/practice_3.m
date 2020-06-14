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

% plot figure

% plot h
% figure();
% histogram(q_h);
% title("param quantize");

% plot distribution

figure();title("non qunatize vs quantize");
histogram(filt_sig);
hold on;
histogram(q_direct_filt_sig);
hold on;
histogram(q_tran_filt_sig);
hold on;
histogram(q_hybrid_filt_sig);
legend("non quantize"," direct quantize", " trans quantize", "hybrid quantize")


sqnr_direct = SQNR(filt_sig, q_direct_filt_sig)
sqnr_tran = SQNR(filt_sig, q_tran_filt_sig)
sqnr_hybrid = SQNR(filt_sig, q_hybrid_filt_sig)


srrc = srrc_pulse(2,2,0.3);
srrc_delay = (length(srrc)-1)/2;

% srrc
srrc_filt_sig = conv(q_val, srrc);
srrc_filt_sig = srrc_filt_sig(srrc_delay+1:end-srrc_delay-1);

srrc_direct_filt = srrc_direct_fix_operation(q_val,srrc);
srrc_tran_filt = srrc_tran_fix_operation(q_val,srrc);
srrc_hybrid_filt = srrc_hybrid_fix_operation(q_val,srrc);

figure();title("srrc filt")
% histogram(srrc_filt_sig);
% hold on;
histogram(srrc_direct_filt);
hold on;
histogram(srrc_tran_filt);
hold on;
histogram(srrc_hybrid_filt);
% "pure srrc",
legend("srrc direct", "srrc tran", "srrc hybrid");




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

function Sout = srrc_direct_fix_operation(x,srrc)

  figure();
  subplot(3,6,1);
  histogram(x);
  title("Sin quantize");

  Sin = srrc(1) * x;
  S1 = real_ADC(Sin, 6, [min(Sin),max(Sin)]);
  subplot(3,6,2);
  histogram(S1);
  title("S1 quantize");

  % first stage
  z1 = srrc(2) * [0 x];
  S2 = real_ADC(z1, 6, [min(z1),max(z1)]);
  subplot(3,6,3);
  histogram(S2);
  title("S2 quantize");

  S5 = [S1 0] + S2;
  S5 = real_ADC(S5, 5, [min(S5),max(S5)]);
  subplot(3,6,4);
  histogram(S5);
  title("S5 quantize");

  % second stage
  z2 = srrc(3) * [0 0 x];
  S3 = real_ADC(z2, 6, [min(z2),max(z2)]);
  subplot(3,6,5);
  histogram(S3);
  title("S3 quantize");

  S6 = [S5 0] + S3;
  S6 = real_ADC(S6, 5, [min(S6),max(S6)]);
  subplot(3,6,6);
  histogram(S6);
  title("S6 quantize");

  % third stage
  z3 = srrc(4) * [0 0 0 x];
  S4 = real_ADC(z3, 6, [min(z3),max(z3)]);
  subplot(3,6,7);
  histogram(S4);
  title("S4 quantize");

  S7 = [S6 0] + S4;
  S7 = real_ADC(S7, 5, [min(S7),max(S7)]);
  subplot(3,6,8);
  histogram(S7);
  title("S7 Direct");

  % 4th stage
  z4 = srrc(5) * [0 0 0 0 x];
  S5 = real_ADC(z4, 6, [min(z4),max(z4)]);
  subplot(3,6,9);
  histogram(S5);
  title("S5 quantize");

  S8 = [S7 0] + S5;
  S8 = real_ADC(S8, 5, [min(S8),max(S8)]);
  subplot(3,6,10);
  histogram(S8);
  title("S8 Direct");

  % 5th stage
  z5 = srrc(6) * [0 0 0 0 0 x];
  S6 = real_ADC(z5, 6, [min(z5),max(z5)]);
  subplot(3,6,11);
  histogram(S6);
  title("S6 quantize");

  S9 = [S8 0] + S6;
  S9 = real_ADC(S9, 5, [min(S9),max(S9)]);
  subplot(3,6,12);
  histogram(S9);
  title("S9 Direct");

  % 6th stage
  z6 = srrc(7) * [0 0 0 0 0 0 x];
  S7 = real_ADC(z6, 6, [min(z6),max(z6)]);
  subplot(3,6,13);
  histogram(S7);
  title("S7 quantize");

  S10 = [S9 0] + S7;
  S10 = real_ADC(S10, 5, [min(S10),max(S10)]);
  subplot(3,6,14);
  histogram(S10);
  title("S10 Direct");

  % 7th stage
  z7 = srrc(8) * [0 0 0 0 0 0 0 x];
  S8 = real_ADC(z7, 6, [min(z7),max(z7)]);
  subplot(3,6, 15);
  histogram(S8);
  title("S8 quantize");

  S11 = [S10 0] + S8;
  S11 = real_ADC(S11, 5, [min(S11),max(S11)]);
  subplot(3,6,16);
  histogram(S11);
  title("S11 Direct");

  % 8th stage
  z8 = srrc(9) * [0 0 0 0 0 0 0 0 x];
  S9 = real_ADC(z8, 6, [min(z8),max(z8)]);
  subplot(3,6, 17);
  histogram(S9);
  title("S9 quantize");

  Sout = [S11 0] + S9;
  Sout = real_ADC(Sout, 5, [min(Sout),max(Sout)]);
  subplot(3,6,18);
  histogram(Sout);
  title("Sout srrc Direct");

end


function Sout = srrc_tran_fix_operation(x,srrc)

  figure();
  subplot(3,6,1);
  histogram(x);
  title("Sin quantize");

  % first stage
  S1 = srrc(1) * x;
  S1 = real_ADC(S1, 6, [min(S1),max(S1)]);
  subplot(3,6,2);
  histogram(S1);
  title("S1 quantize");

  z1 = [0 S1];
  S2 = real_ADC(z1, 5, [min(z1),max(z1)]);
  subplot(3,6,3);
  histogram(S2);
  title("S2 quantize");


  S3 = srrc(2)*x;
  S3 = real_ADC(S3, 5, [min(S3),max(S3)]);
  subplot(3,6,4);
  histogram(S3);
  title("S3 quantize");

  S4 = z1+[S3 0];
  S4 = real_ADC(S4, 5, [min(S4),max(S4)]);
  subplot(3,6,5);
  histogram(S4);
  title("S4 quantize");

  % second stage
  S5 = srrc(3)*x;
  S5 = real_ADC(S5, 5, [min(S5),max(S5)]);
  subplot(3,6,6);
  histogram(S5);
  title("S5 quantize");

  z2 = [0 S4];
  S6 = z2+[S5 0 0];
  subplot(3,6,7);
  histogram(S6);
  title("S6 quantize");

  % third stage
  S7 = srrc(4)*x;
  S7 = real_ADC(S7, 5, [min(S7),max(S7)]);
  subplot(3,6,7);
  histogram(S7);
  title("S7 quantize");

  z3 = [0 S6];
  S8 = z3+[S7 0 0 0];
  subplot(3,6,8);
  histogram(S8);
  title("S8 transpose");

  % 4th stage
  S9 = srrc(5)*x;
  S9 = real_ADC(S9, 5, [min(S9),max(S9)]);
  subplot(3,6,7);
  histogram(S9);
  title("S9 quantize");

  z4 = [0 S8];
  S10 = z4+[S9 0 0 0 0];
  subplot(3,6,8);
  histogram(S10);
  title("S10 transpose");

  % 5th stage
  S11 = srrc(5)*x;
  S11 = real_ADC(S11, 5, [min(S11),max(S11)]);
  subplot(3,6,9);
  histogram(S11);
  title("S11 quantize");

  z5 = [0 S10];
  S12 = z5+[S11 0 0 0 0 0];
  subplot(3,6,10);
  histogram(S12);
  title("S12 transpose");

  % 6th stage
  S13 = srrc(6)*x;
  S13 = real_ADC(S13, 5, [min(S13),max(S13)]);
  subplot(3,6,9);
  histogram(S13);
  title("S13 quantize");

  z6 = [0 S12];
  S14 = z6+[S13 0 0 0 0 0 0];
  subplot(3,6,10);
  histogram(S14);
  title("S14 transpose");

  % 7th stage
  S15 = srrc(7)*x;
  S15 = real_ADC(S15, 5, [min(S15),max(S15)]);
  subplot(3,6,11);
  histogram(S15);
  title("S15 quantize");

  z7 = [0 S14];
  S16 = z7+[S15 0 0 0 0 0 0 0];
  subplot(3,6,12);
  histogram(S16);
  title("S16 transpose");

  % 8th stage
  S17 = srrc(8)*x;
  S17 = real_ADC(S17, 5, [min(S17),max(S17)]);
  subplot(3,6,13);
  histogram(S17);
  title("S17 quantize");

  z8 = [0 S16];
  S18 = z8+ [S17 0 0 0 0 0 0 0 0];
  subplot(3,6,14);
  histogram(S18);
  title("S18 transpose");

  % 9th stage
  S19 = srrc(9)*x;
  S19 = real_ADC(S19, 5, [min(S19),max(S19)]);
  subplot(3,6,15);
  histogram(S19);
  title("S19 quantize");

  z9 = [0 S18];
  Sout = z9+ [S19 0 0 0 0 0 0 0 0 0];
  subplot(3,6,16);
  histogram(Sout);
  title("Sout tran transpose");
  srrc_delay = (length(srrc)-1)/2;
  Sout = conv(x,srrc);
  Sout = Sout(srrc_delay+1:end-srrc_delay);
end

function Sout = srrc_hybrid_fix_operation(x,srrc)

  figure();
  subplot(3,6,1);
  histogram(x);
  title("Sin quantize");
  z1 = [0 x];
  % zn = z1;
  % ori_sig = srrc(1)*[x 0 0]+srrc(2)*[z1 0]+srrc(3)*[0 z1];
  % for i = 1:6
  %   zn = [0 zn];
  %   first_two_stage(ori_sig, srrc);
  % end
  % first stage
  S1 = srrc(1) * x;
  S1 = real_ADC(S1, 6, [min(S1),max(S1)]);
  subplot(3,6,2);
  histogram(S1);
  title("S1 quantize");

  z1 = [0 S1];
  S2 = real_ADC(z1, 5, [min(z1),max(z1)]);
  subplot(3,6,3);
  histogram(S2);
  title("S2 quantize");


  S3 = srrc(2)*x;
  S3 = real_ADC(S3, 5, [min(S3),max(S3)]);
  subplot(3,6,4);
  histogram(S3);
  title("S3 quantize");

  S4 = z1+[S3 0];
  S4 = real_ADC(S4, 5, [min(S4),max(S4)]);
  subplot(3,6,5);
  histogram(S4);
  title("S4 quantize");

  % second stage
  S5 = srrc(3)*x;
  S5 = real_ADC(S5, 5, [min(S5),max(S5)]);
  subplot(3,6,6);
  histogram(S5);
  title("S5 quantize");

  z2 = [0 S4];
  S6 = z2+[S5 0 0];
  subplot(3,6,7);
  histogram(S6);
  title("S6 quantize");

  % third stage
  S7 = srrc(4)*x;
  S7 = real_ADC(S7, 5, [min(S7),max(S7)]);
  subplot(3,6,7);
  histogram(S7);
  title("S7 quantize");

  z3 = [0 S6];
  S8 = z3+[S7 0 0 0];
  subplot(3,6,8);
  histogram(S8);
  title("S8 transpose");

  % 4th stage
  S9 = srrc(5)*x;
  S9 = real_ADC(S9, 5, [min(S9),max(S9)]);
  subplot(3,6,7);
  histogram(S9);
  title("S9 quantize");

  z4 = [0 S8];
  S10 = z4+[S9 0 0 0 0];
  subplot(3,6,8);
  histogram(S10);
  title("S10 transpose");

  % 5th stage
  S11 = srrc(5)*x;
  S11 = real_ADC(S11, 5, [min(S11),max(S11)]);
  subplot(3,6,9);
  histogram(S11);
  title("S11 quantize");

  z5 = [0 S10];
  S12 = z5+[S11 0 0 0 0 0];
  subplot(3,6,10);
  histogram(S12);
  title("S12 transpose");

  % 6th stage
  S13 = srrc(6)*x;
  S13 = real_ADC(S13, 5, [min(S13),max(S13)]);
  subplot(3,6,9);
  histogram(S13);
  title("S13 quantize");

  z6 = [0 S12];
  S14 = z6+[S13 0 0 0 0 0 0];
  subplot(3,6,10);
  histogram(S14);
  title("S14 transpose");

  % 7th stage
  S15 = srrc(7)*x;
  S15 = real_ADC(S15, 5, [min(S15),max(S15)]);
  subplot(3,6,11);
  histogram(S15);
  title("S15 quantize");

  z7 = [0 S14];
  S16 = z7+[S15 0 0 0 0 0 0 0];
  subplot(3,6,12);
  histogram(S16);
  title("S16 transpose");

  % 8th stage
  S17 = srrc(8)*x;
  S17 = real_ADC(S17, 5, [min(S17),max(S17)]);
  subplot(3,6,13);
  histogram(S17);
  title("S17 quantize");

  z8 = [0 S16];
  S18 = z8+ [S17 0 0 0 0 0 0 0 0];
  subplot(3,6,14);
  histogram(S18);
  title("S18 transpose");

  % 9th stage
  S19 = srrc(9)*x;
  S19 = real_ADC(S19, 5, [min(S19),max(S19)]);
  subplot(3,6,15);
  histogram(S19);
  title("S19 quantize");

  z9 = [0 S18];
  Sout = z9+ [S19 0 0 0 0 0 0 0 0 0];
  subplot(3,6,16);
  histogram(Sout);
  title("Sout tran transpose");
  srrc_delay = (length(srrc)-1)/2;
  Sout = conv(x,srrc);
  Sout = Sout(srrc_delay+1:end-srrc_delay);
end


function out = first_two_stage(ori_sig,sig,srrc)
  % zn = [0 ori_sig];
  conpensate = length(sig)-length(ori_sig);
  out = [ori_sig zeros(1,conpensate)] + sig;
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
