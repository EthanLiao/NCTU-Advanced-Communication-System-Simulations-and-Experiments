% Parametors
clf;close all;close all;
fs = 1*10^6;
fd = 150*10^3;
fIF = 2*10^6;
fc = 8*10^6;
BT = 0.5;
f_DAC = 16*10^6;
f_DMA = 64*10^6;
signal = [-1,1,-1,1];
L = 16;

gau_filter = Gfilter(BT,f_DAC/fs);
gau_delay = (length(gau_filter)-1)/2;

srrc_16 = srrc_pulse(f_DAC/fs,5,0.3);
srrc_16_delay = (length(srrc_16)-1)/2;

srrc_64 = srrc_pulse(f_DMA/f_DAC,5,0.3);
srrc_64_delay = (length(srrc_64)-1)/2;

% Gaussian filter
t_DAC_sig = DAC(signal, f_DAC/fs);
g_sig = conv(t_DAC_sig, gau_filter);
g_sig = g_sig(gau_delay+1:end-gau_delay);

% sum
for i = 1:length(g_sig)
  sum_gsig(i) = sum(g_sig(1:i));
end
cp_sig = exp(1j*2*pi*fd/fs*1/fs*sum_gsig);



% IF modulation
t = [0:length(cp_sig)-1];
IF_sig = real(cp_sig.*exp(1j*2*pi*fIF/fs*t));


% DAC/F
t_DMA_sig = conv(DAC(IF_sig,f_DMA/f_DAC),srrc_64);
t_DMA_sig = t_DMA_sig(srrc_64_delay+1:end-srrc_64_delay);

% F/ADC
r_DMA_sig = conv(t_DMA_sig,srrc_64);
r_DMA_sig = r_DMA_sig(srrc_64_delay+1:end-srrc_64_delay);
r_DMA_sig = ADC(r_DMA_sig, f_DMA/f_DAC);

% IF demodulation
t = [0:length(r_DMA_sig)-1];
dmod_IF_sig = r_DMA_sig .* exp(-1j*2*pi*fIF/fs*t);

% F/PH
re_srrc = conv(dmod_IF_sig,srrc_16);
re_srrc = re_srrc(srrc_16_delay+1:end-srrc_16_delay);

% phase and diff
re_phase = phase(re_srrc);
re_phase_diff(1) = re_phase(1)/(2*pi*fd/fs*1/fs);
for i = 2:length(re_phase)
  re_phase_diff(i) = (re_phase(i)-re_phase(i-1))/(2*pi*fd/fs*1/fs);
end

% Gaussian Filter and ADC
r_sig = conv(re_phase_diff, gau_filter);
r_sig = r_sig(gau_delay+1:end-gau_delay);
r_sig = ADC(r_sig, f_DAC/fs);

r_sig = r_sig /max(abs(r_sig));
figure()
stem(r_sig);title('Recieved signal');

% some function for reuse
function up_sig = DAC(sig,up_factor)
  up_sig = zeros(1,length(sig)*up_factor);
  up_sig(1:up_factor:end) = sig;
end

function down_sig = ADC(sig, down_factor)
  down_sig = zeros(1, length(sig));
  down_sig = sig(1:down_factor:end);
end

% Gaussian Filter
function g_filter = Gfilter(BT,M)
  t = [-16:16];
  g_filter = sqrt(2*pi/log(2))*exp(-2*(pi^2)/log(2)*(BT/M)^2*(t.^2));
ends


function [y,t] = srrc_pulse(T,A,a)
  t = [-A*T:A*T] + 10^(-8);
  if (a>0 && a<=1)
    num = cos((1+a)*pi*t/T) + T*sin((1-a)*pi*t/T)./(4*a*t);
    denom = 1-(4*a*t/T).^2;
    y = (4*a/pi) * num./denom;
  else
    y = 1/T * sin(pi.*t./T) ./ (pi*t./T);
  end
end


% % carrier modulation
% t = [0:length(t_DMA_sig)-1];
% mod_sig = DMA_sig.*exp(1j*2*pi*(fc-fIF)/fs*t);
%
%
% % carrier demodulation
% t = [0:length(mod_sig)-1];
% dmod_sig = mod_sig.*exp(-1j*2*pi*(fc-fIF)/fs*t);
