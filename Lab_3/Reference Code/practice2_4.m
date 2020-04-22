clf
N = 512
x = linspace(-N/2,N/2)
y= sinc(x)
plot(x,abs(fft(y)))
