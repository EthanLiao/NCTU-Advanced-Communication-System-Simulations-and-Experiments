clf
% ------rectangular wave------
ele_nums = 1024
mid = ceil(ele_nums/2)
half = 5

% rectangular pulse
x = zeros(1, ele_nums)
x(mid-half : mid+half) = 1
subplot(3,1,1);stem(abs(fft(rect)),'.-');title('rectangular wave freq response');grid on;
% ------filter the wave------


% -----2-2HP----
C = [1, 0.9999]
% 補償項
y = zeros(1, ele_nums+length(C)-1)

% 做convolution的項次
for n = 2:ele_nums
  y(n) = x(n)-C(2)*y(n-1)
end
subplot(3,1,2);stem(abs(fft(y)),'.-');title('HP frequency response');grid on;

% -----2-2LP----
L_C = [1, 0.16, 0.65]
% 補償項
L_y = zeros(1, ele_nums+length(L_C)-1)
% 做convolution的項次
for n = 3:ele_nums
  L_y(n) = x(n)+L_C(2)*L_y(n-1)-L_C(3)*L_y(n-2)
end
subplot(3,1,3);stem(abs(fft(L_y)),'.-');title('LP frequency response');grid on;
