clf;clear all;close all;
BT = 0.5;
pulse = Gfilter(BT, 10);
subplot(2,1,1);plot(pulse);title("time domain Gaussian filter");
subplot(2,1,2);plot(abs(fftshift(pulse)));title("frequncy domain Gaussian filter");

function g_filter = Gfilter(B,L)
  t = [-L:0.001:L];
  g_filter = sqrt(2*pi/log10(2))*B*exp(-2*pi/log10(2)*B^2*(t.^2));
end
