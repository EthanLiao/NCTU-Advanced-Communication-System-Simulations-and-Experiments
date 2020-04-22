clf
% ------rectangular wave------
ele_nums = 16
mid = ceil(ele_nums/2)
half = 5

% rectangular pulse
x = zeros(1, ele_nums)
x(mid-half : mid+half) = 1
subplot(3,1,1);plot(abs(fft(rect)),'.-');title('rectangular wave freq response');grid on;
% ------filter the wave------


% -----2-1HP----
C = [1, -0.9]
% 補償項
y(1) = C(1) * x(1)

% 做convolution的項次
for n = 2:ele_nums
  y(n) = C(1) * x(n) + C(2)*x(n-1)
end
subplot(3,1,2);plot(abs(fft(y)),'.-');title('HP_freq_response');grid on;

% % -----2-1LP----
C = [1, 0.9]
% 補償項
hp_y(1) = C(1) * x(1)

% 做convolution的項次
for n = 2:ele_nums
  hp_y(n) = C(1) * x(n) + C(2)*x(n-1)
end
subplot(3,1,3);plot(abs(fft(hp_y)),'.-');title('LP_freq_response');grid on;
