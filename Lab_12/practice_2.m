clf;clear all;close all;
BT = 0.5;
pulse = Gfilter(BT, 10);
subplot(2,1,1);plot(pulse);title("time domain Gaussian filter");
subplot(2,1,2);plot(abs(fftshift(pulse)));title("frequncy domain Gaussian filter");

function g_filter = Gfilter(BT,M)
  t = [-16:16];
  g_filter = sqrt(2*pi/log(2))*exp(-2*(pi^2)/log(2)*(BT/M)^2*(t.^2));
end
