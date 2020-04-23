clf;clear all;
srrc = srrc_pulse(5,10,0.3)
rc = conv(srrc,srrc)


subplot(2,2,1);plot(srrc);title('srrc sginal')
subplot(2,2,2);plot(fftshift(abs(fft(srrc))));title('srrc sginal f domain')
subplot(2,2,3);plot(rc);title('rc sginal')
subplot(2,2,4);plot(fftshift(abs(fft(rc))));title('rc sginal f domain')

function [y,t] = srrc_pulse(T,A,a)
  t = [-A*T:A*T]+10^(-8)
  if (a>0 && a<1)
    num = cos((1+a)*pi*t./T)+T*sin((1-a)*pi*t./T)./(4*a*t)
    denom = 1-(4*a*t./T).^2
    y = 4*a/pi.*num./denom
  else
    y = 1/T * sin(pi*t./T)./(pi*t./T)
  end
end
