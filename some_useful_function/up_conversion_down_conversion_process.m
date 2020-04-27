rec_len = 100
mid = ceil(rec_len/2)
half = 20
rect = zeros(1,rec_len)
rect(mid-half:mid+half) = 1
tri_pulse = conv(rect,rect)
signal = [-1,1,-1,1,-1,-1,1]
srrc_4 = srrc_pulse(4,10,0.3)
DAC_sampling_factor = 4*10^6 / 10^6
DMA_sampling_factor = 32*10^6 / (4*10^6)
carrier_freq = 8*10^8
load('./filter/conversion_filter')

trans_sig = filter(conversion_filter,DAC(conv(DAC(signal,DAC_sampling_factor),srrc_4,'same'),DMA_sampling_factor))

t = 1:length(trans_sig)
carrier = cos(2*pi*carrier_freq*t)+i*sin(2*pi*carrier_freq*t)
mod_sig = real(trans_sig.*carrier)
dmod_sig = mod_sig.*conj(carrier)

rcv_sig = real(ADC(conv(ADC(filter(conversion_filter,dmod_sig),DMA_sampling_factor),srrc_4,'same'),DAC_sampling_factor))

subplot(2,1,1);stem(signal)
subplot(2,1,2);stem(rcv_sig./max(rcv_sig))

function srrc = srrc_pulse(T,A,a)
  t = [-A*T:A*T]+10^(-8)
  if(a>0 && a<=1)
    num = cos((1+a)*pi*t/T)+sin((1-a)*pi*t/T)./(4*a*t/T)
    denum = 1-(4*a*t/T).^2
    srrc = 4*a/pi*num./denum
  else
    srrc = 1/T * sin(pi*t/T)./(pi*t/T)
  end
end

function down_sig = ADC(sig,down_factor)
  down_sig = zeros(1,length(sig))
  down_sig = sig(1:down_factor:end)
end

function up_sig = DAC(sig,up_factor)
  up_sig = zeros(1,length(sig)*up_factor)
  up_sig(1:up_factor:end) = sig
end

function y = add_awgn_noise(sig,SNR_DB)
  L = length(sig)
  snr = 10^(SNR_DB/10)
  SYME = sum(abs(sig).^2)/L
  N0 = SYME / snr
  if isreal(x)
    n = sqrt(N0)*randn(1,L)
  else
    n = N0/2*(randn(1,L)+i*randn(1,L))
  y = x+n
  end
end
