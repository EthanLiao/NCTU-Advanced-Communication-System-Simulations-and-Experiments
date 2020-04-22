clf
N = 256
fc_q = 1/4
fc_oc = 1/8
time = (0:N-1)
x = cos(2*pi*fc_q*time) + cos(2*pi*fc_oc*time)
% C = [1, -1.848, 1.000465]
% 補�?��??

% ------rectangular wave------


% y = zeros(1, N+length(C)-1)
% y(1) = a(1) * x(1)
% y(2) = a(1) * x(2) + a(2)*x(1)
%
% % ??�convolution??��?�次
% for n = 3:N
%   y(n) = a(1) * x(n) + a(2)*x(n-1) + a(3)*x(n-2) - b(2)*y(n-1)- b(3)*y(n-2)
% end
% -------------
a = poly([0.95999*(cos(1/4*pi)+i*sin(1/4*pi)), 0.95999*(conj(cos(1/4*pi)+i*sin(1/4*pi)))])
b = poly([1.09656*(cos(1/8*pi)+i*sin(1/8*pi)), 1.09656*(conj(cos(1/8*pi)+i*sin(1/8*pi)))])

[h,w] = freqz([b,a,1024]);
subplot(4,1,1);plot(w, 360/(2*pi)*angle(h))
SlopeBetweenPoints = diff(360/(2*pi)*angle(h))./diff(w)
SlopeBetweenPoints = reshape(SlopeBetweenPoints,1,length(SlopeBetweenPoints))
delay = SlopeBetweenPoints(1)
filter_sig = 0.12*filter(b,a,x)
filter_sig = filter_sig(150:202)
sig = cos(2*pi*fc_oc*time+delay)
sig = sig(150:202)
% indicies = 0:1:-delay
% test(indicies)=[]
subplot(4,1,2);plot(filter_sig)
subplot(4,1,3);plot(sig)
IIR_snr = 10*log10(moment(sig,2) / moment(abs(filter_sig-sig),2))
% --------------
% mid = ceil(N/2)
% half = 100
% rect = zeros(1, N)
% rect(mid-half : mid+half) = 1
% freq_y = zeros(1, N+length(C)-1)
% % freq_y(1) = C(1) * rect(1)
% % freq_y(2) = C(1) * rect(2) + C(2)*rect(1)
%
% % % ??�convolution??��?�次
% for n = 3:N
%   freq_y(n) = C(1) * rect(n) + C(2)*rect(n-1) + C(3)*rect(n-2) - 0.36*freq_y(n-2)
% end


% subplot(4,1,1);plot(time,x,'.-');title('sinusoid signal combined');grid on;
% subplot(4,1,2);plot(y,'.-');title('sinusoid signal filtered');grid on;
% subplot(4,1,3);plot(cos(2*pi*fc_q*time),'.-');title('origin');grid on;
% subplot(4,1,4);plot(abs(fft(freq_y)),'.-');title('freq response');grid on;
