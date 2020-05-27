clf;clear all;close all;
BT = 0.5;
M = 17;
fb = 2;
pulse = Gfilter(BT, M, fb);
subplot(2,1,1);plot(pulse);title("time domain Gaussian filter");
subplot(2,1,2);plot(abs(fftshift(pulse)));title("frequncy domain Gaussian filter");
BT = [0.2, 0.25, 0.3];
for i = 1:length(BT)
  plot(Gfilter(BT(i),M,1/fb));
  hold on;
end
legend("WTb = 0.2","WTb = 0.25", "WTb = 0.3" )
% Gaussian Filter
function g_filter = Gfilter(BT, M, Tb)
  t = [-64:64];
  B = BT/Tb;
  C = sqrt(2*pi/log(2)) * B;
  g_filter = C*exp(-2*(pi^2)/log(2)*(BT/M)^2*(t.^2));
  g_filter = g_filter ./ sqrt(sum(g_filter.^2));
  % g_filter = g_filter ./ max(g_filter.^2);
end
