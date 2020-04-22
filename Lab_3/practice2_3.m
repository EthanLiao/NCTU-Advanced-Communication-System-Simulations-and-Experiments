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


% -----2-2BP----
X_C = [1, -1.6, -0.16]
Y_C = [0, 0, 0.81]
% 補償項
y = zeros(1, ele_nums+length(C)-1)

% 做convolution的項次
for n = 3:ele_nums
  y(n) = X_C(1)*x(n)+X_C(2)*x(n-2)-Y_C(3)*y(n-2)
end
subplot(3,1,2);stem(abs(fft(y)),'.-');title('BP_freq_response');grid on;
