clf
N = 256
fc_q = 1/4
fc_oc = 1/8
time = (0:N-1)
x = cos(2*pi*fc_q*time) + cos(2*pi*fc_oc*time)
pole = poly([0.95999*(cos(1/4*pi)+i*sin(1/4*pi)), 0.95999*(conj(cos(1/4*pi)+i*sin(1/4*pi)))])
zero = poly([1.09656*(cos(1/8*pi)+i*sin(1/8*pi)), 1.09656*(conj(cos(1/8*pi)+i*sin(1/8*pi)))])

[h,w] = freqz([zero,pole,1024]);
subplot(4,1,1);plot(w, 360/(2*pi)*angle(h))
SlopeBetweenPoints = diff(360/(2*pi)*angle(h))./diff(w)
SlopeBetweenPoints = reshape(SlopeBetweenPoints,1,length(SlopeBetweenPoints))
delay = SlopeBetweenPoints(1)
filter_sig = 0.12*filter(zero,pole,x)
filter_sig = filter_sig(150:202)
sig = cos(2*pi*fc_oc*time+delay)
sig = sig(150:202)


subplot(4,1,2);plot(filter_sig)
subplot(4,1,3);plot(sig)
IIR_snr = 10*log10(moment(sig,2) / moment(abs(filter_sig-sig),2))
