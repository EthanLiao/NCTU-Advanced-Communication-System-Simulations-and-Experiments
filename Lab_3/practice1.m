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
% subplot(3,1,1);stem(rect,'.-');title('rectangular wave');grid on;
% subplot(3,1,2);stem(y,'.-');title('response');grid on;
% subplot(3,1,3);plot(abs(fft(y)),'.-');title('response');grid on;

x = [1 2 3 ]
h = [4 3 2]
xs = length(x)
hs = length(h)
y = zeros(1,xs+hs-1)
xe = [zeros(1,hs-1) x zeros(1,hs-1)]
for k =1:xs-hs+1
  y(k) = h(end:-1:1)*x(k:k+hs-1).'
end
stem(y)
