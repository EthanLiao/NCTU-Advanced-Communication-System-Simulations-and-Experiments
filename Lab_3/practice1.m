clf;clear all;
% % ------rectangular wave------
% ele_nums = 1024
% mid = ceil(ele_nums/2)
% half = 10
%
% % rectangular pulse
% x = zeros(1, ele_nums)
% x(mid-half : mid+half) = 1
%
% % ------filter the wave------
%
% C = [1, -1.848, 1.000465]
% % 補償項
% y(1) = C(1) * x(1)
% y(2) = C(1) * x(2) + C(2)*x(1)
%
% % 做convolution的項次
% for n = 3:ele_nums
%   y(n) = C(1) * x(n) + C(2)*x(n-1) + C(3)*x(n-2)
% end
%
%
% % ------plot the figure------
% figure()
% subplot(3,1,1);stem(rect,'.-');title('rectangular wave');grid on;
% subplot(3,1,2);stem(y,'.-');title('response');grid on;
% subplot(3,1,3);plot(abs(fft(y)),'.-');title('response');grid on;

% practice 1 teacher's version
% x = [1 2 3 4]
% h = [4 3 2]
% xs = length(x)
% hs = length(h)
% y = zeros(1,xs+hs-1)
% xe = [zeros(1,hs-1) x zeros(1,hs-1)]
% for k =1:xs-hs+1
%   y(k) = h(end:-1:1)*x(k:k+hs-1).'
% end
% stem(y)

% practice1 _version 2

fc_6 = 1/6
fc_8 = 1/8
N = 10
t_6 = (0:0.01:N)/fc_6
t_8 = (0:0.01:N)/fc_8
signal_6 = sin(2*pi*fc_6*t_6)
signal_8 = cos(2*pi*fc_8*t_8)
% signal_8 = zeros(1,length(signal_6))
% signal_8(1:end) = 1
combined_signal = signal_6+signal_8

% FIR filter
zero = poly(0.9959*[cos(pi*fc_6)+i*sin(pi*fc_6),cos(pi*fc_6)-i*sin(pi*fc_6)])
pole = 1
FIR_signal = filter(zero,pole,combined_signal)
% IIR filter
zero = poly([0.9959*(cos(pi*fc_6)+j*sin(pi*fc_6)),0.9959*(cos(pi*fc_6)-j*sin(pi*fc_6))])
pole = poly([0.9*(cos(pi*fc_8)+j*sin(pi*fc_8)),0.9*(cos(pi*fc_8)-j*sin(pi*fc_8))])
IIR_signal = filter(zero,pole,combined_signal)


subplot(5,1,1);plot(delay(signal_8,31));title('original signal')
subplot(5,1,2);plot(signal_6);title('another signal')
subplot(5,1,3);plot(combined_signal);title('combined signal')
subplot(5,1,4);plot(FIR_signal./0.37169);title('FIR filtered signal')
subplot(5,1,5);plot(IIR_signal./max(IIR_signal));title('IIR filtered signal')

function delayed_signal = delay(sig,delay_num)
  delta = zeros(1,2*delay_num)
  delta(delay_num) = 1
  delayed_signal = conv(sig,delta)
end
