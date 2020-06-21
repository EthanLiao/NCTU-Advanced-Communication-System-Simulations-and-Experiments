close all;clf;clear all;

%%%%%% question 2 and 4 %%%%%%
N = 2048;
sig = randi([0 1],1,N);
sig(sig==0) = -1;

c_sig = transbranch(sig);
r_DAC_sig = rcvbranch(c_sig);

snr_db = SNR(sig, r_DAC_sig)

subplot(2,1,1);stem(sig);title("transmitted signal");
subplot(2,1,2);stem(r_DAC_sig);title("recieved signal");

%%%%%% question 5 %%%%%%
N = 20
tap = zeros(1,19);
signal = [];
for i = 1:N
  signal = [signal 1 tap];
end

multipath = [1 tap -1.2 tap -0.25 tap 0.3];


c_sig = ISI_transbranch(signal);
c_sig = conv(c_sig, multipath);
r_DAC_sig = ISI_rcvbranch(c_sig);

equalizer = 20*(-1/56)*((1/2).^[0:length(r_DAC_sig)-1]);
equalized_sig = conv(r_DAC_sig, equalizer);

equalizer = 20*(36/119/5)*6/5*((6/5).^[-length(equalized_sig):-1]);
equalized_sig = conv(equalized_sig, equalizer);

equalizer = 20*(1/68/2)*((-1/2).^[0:length(equalized_sig)-1]);
equalized_sig = conv(equalized_sig, equalizer);
equalized_sig = equalized_sig / max(equalized_sig);

group_delay = 806;
equalized_sig = equalized_sig(group_delay:end);
equalized_sig = equalized_sig(1:400);
snr_db_equ = SNR(signal, equalized_sig)

figure()
subplot(2,1,1);stem(signal);title("transmitted signal");
subplot(2,1,2);stem(equalized_sig); title("equalized signal");

function c_sig = transbranch(sig)
  % Some parametors
  fb = 1e6;
  freq_DAC = 16e6;
  freq_DMA = 128e6;
  fc = 32e6;
  f_DAC = freq_DAC / fb;
  f_DMA = freq_DMA / freq_DAC;
  load('./filter/IIR_filter');
  group_delay = 7;
  srrc_16 = rcosine(1, 16, 'fir/sqrt', 0.3, 5);
  srrc_delay = (length(srrc_16)-1)/2;

  % DMA
  t_DAC_sig = conv(srrc_16, DAC(sig, f_DAC));
  t_DAC_sig = t_DAC_sig(srrc_delay+1:end-srrc_delay);
  t_DAC_sig = t_DAC_sig / max(t_DAC_sig);

  % DAC
  t_DMA_sig = filter(IIR_filter, DAC(t_DAC_sig, f_DMA));
  t_DMA_sig = t_DMA_sig(group_delay:end);
  t_DMA_sig = t_DMA_sig / max(t_DMA_sig);

  c_sig = real(t_DMA_sig .* exp(1j*pi*fc/freq_DMA*[0:length(t_DMA_sig)-1]));
end

function r_DAC_sig = rcvbranch(c_sig)
  % Some parametors
  fb = 1e6;
  freq_DAC = 16e6;
  freq_DMA = 128e6;
  fIF = 2e6;
  fc = 32e6;
  f_DAC = freq_DAC / fb;
  f_DMA = freq_DMA / freq_DAC;
  load('./filter/IIR_filter');
  group_delay = 7;
  srrc_16 = rcosine(1, 16, 'fir/sqrt', 0.3, 5);
  srrc_delay = (length(srrc_16)-1)/2;

  % transfer to the IF band
  c_sig = c_sig .* cos(2*pi*(fc-fIF)/freq_DMA*[0:length(c_sig)-1]);

  % DMA
  r_DMA_sig = filter(IIR_filter, c_sig);
  r_DMA_sig = r_DMA_sig(group_delay:end);
  r_DMA_sig = ADC(r_DMA_sig, f_DMA);
  r_DMA_sig = r_DMA_sig / max(r_DMA_sig);

  % transfer to the baseband
  r_DMA_sig = real(r_DMA_sig .* exp(-1j*2*pi*fIF/freq_DAC*[0:length(r_DMA_sig)-1]));

  % DAC
  r_DAC_sig = conv(srrc_16, r_DMA_sig);
  r_DAC_sig = r_DAC_sig(srrc_delay+1:end-srrc_delay);
  r_DAC_sig = ADC(r_DAC_sig, f_DAC);
  r_DAC_sig = r_DAC_sig / max(r_DAC_sig);
end


function c_sig = ISI_transbranch(sig)
  % Some parametors
  freq_DMA = 128e6;
  fc = 32e6;
  f_DAC = 3;
  f_DMA = 8;
  srrc_3 = rcosine(1, 3, 'fir/sqrt', 0.3, 5);
  srrc_3_delay = (length(srrc_3)-1)/2;

  srrc_8 = rcosine(1, 8, 'fir/sqrt', 0.3, 5);
  srrc_8_delay = (length(srrc_8)-1)/2;

  % DAC
  t_DAC_sig = conv(srrc_3, DAC(sig, f_DAC));
  t_DAC_sig = t_DAC_sig(srrc_3_delay+1:end-srrc_3_delay);
  % t_DAC_sig = t_DAC_sig / max(t_DAC_sig);

  % DMA
  t_DMA_sig = conv(srrc_8, DAC(t_DAC_sig, f_DMA));
  t_DMA_sig = t_DMA_sig(srrc_8_delay+1:end-srrc_8_delay);
  % t_DMA_sig = t_DMA_sig / max(t_DMA_sig);

  c_sig = real(t_DMA_sig .* exp(1j*2*pi*fc/freq_DMA*[0:length(t_DMA_sig)-1]));
end

function r_DAC_sig = ISI_rcvbranch(c_sig)
  % Some parametors
  freq_DMA = 128e6;
  fc = 32e6;
  f_DAC = 3;
  f_DMA = 8;
  srrc_3 = rcosine(1, 3, 'fir/sqrt', 0.3, 5);
  srrc_3_delay = (length(srrc_3)-1)/2;

  srrc_8 = rcosine(1, 8, 'fir/sqrt', 0.3, 5);
  srrc_8_delay = (length(srrc_8)-1)/2;

  c_sig = c_sig .* exp(-1j*2*pi*fc/freq_DMA*[0:length(c_sig)-1]);

  % DMA
  r_DMA_sig = conv(srrc_8,c_sig);
  r_DMA_sig = r_DMA_sig(srrc_8_delay+1:end-srrc_8_delay);
  r_DMA_sig = ADC(r_DMA_sig, f_DMA);
  % r_DMA_sig = r_DMA_sig / max(r_DMA_sig);

  % DAC
  r_DAC_sig = conv(srrc_3,r_DMA_sig);
  r_DAC_sig = r_DAC_sig(srrc_3_delay+1:end-srrc_3_delay);
  r_DAC_sig = ADC(r_DAC_sig, f_DAC);
  % r_DAC_sig = r_DAC_sig / max(r_DAC_sig);

  r_DAC_sig = real(r_DAC_sig);
end


function d_sig = ADC(sig ,d_fac)
  d_sig = zeros(1,length(sig));
  d_sig = sig(1:d_fac:end);
end

function u_sig = DAC(sig, u_fac)
  u_sig = zeros(1,length(sig)*u_fac);
  u_sig(1:u_fac:end) = sig;
end

function snr = SNR(t_sig, r_sig)
  snr = 10*log10(mean(abs(t_sig).^2)/ mean(abs(r_sig-t_sig).^2));
end
