N = 16
fc_q = 1/4
fc_oc = 1/8
time = (0:N-1)
x = cos(2*pi*fc_q*time) + cos(2*pi*fc_oc*time)
C = [1, -1.848, 1.000465]
% 補償項
y(1) = C(1) * x(1)
y(2) = C(1) * x(2) + C(2)*x(1)

% 做convolution的項次
for n = 3:N
  y(n) = C(1) * x(n) + C(2)*x(n-1) + C(3)*x(n-2)
end
subplot(3,1,1);plot(time,x,'.-');title('sinusoid signal combined');grid on;
