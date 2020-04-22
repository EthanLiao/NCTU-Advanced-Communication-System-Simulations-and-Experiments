clf
N = 256
fc_q = 1/4
fc_oc = 1/8
time = (0:N-1)
x = cos(2*pi*fc_q*time) + cos(2*pi*fc_oc*time)
pole = poly([0.7*(cos(1/4*2*pi)+i*sin(1/4*2*pi)), 0.7*(conj(cos(1/4*2*pi)+i*sin(1/4*2*pi)))])
zero = poly([1*(cos(1/8*2*pi)+i*sin(1/8*2*pi)), 1*(conj(cos(1/8*2*pi)+i*sin(1/8*2*pi)))])
[h,w] = freqz([zero,pole,1024]);
SlopeBetweenPoints = diff(360/(2*pi)*angle(h))./diff(w)
SlopeBetweenPoints = reshape(SlopeBetweenPoints,1,length(SlopeBetweenPoints))
delay = SlopeBetweenPoints(1)
filter_sig = 0.12*filter(zero,pole,x)
filter_sig = filter_sig(150:202)
sig = cos(2*pi*fc_q*time+delay)
sig = sig(150:202)
IIR_snr = 10*log10(moment(sig,2) / moment(abs(filter_sig-sig),2))
% zplane(zero,pole)
subplot(4,1,1);plot(x);title('Signal Combined');grid on;
subplot(4,1,2);plot(sig);title('Original Signal');grid on;
subplot(4,1,3);plot(filter_sig);title('IIR fitered signal');grid on;
subplot(4,1,4);plot(w, 360/(2*pi)*angle(h));title('Phase Response');grid on;
% freqz([zero,pole,1024])
