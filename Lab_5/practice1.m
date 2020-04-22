% 8PAM Mapping
M_PAM = 8
x = (0:7)'
y = pammod(x,M_PAM,pi/2)
scatterplot(y)
text(real(y)+0.1,imag(y)-0.1,dec2bin(x))
% 16 QAM Mapping
M_QAM = 16  % modulation order
x = (0:15)' % generate a modulation symbol
y = qammod(x,M_QAM,'gray')
scatterplot(y)
text(real(y)+0.1,imag(y)+0.1,dec2bin(x))
title('16 QAM gray mapping')


N = 3.0;
x=linspace(-N, N);
y=x;
[X,Y]=meshgrid(x,y);
z=(1000/sqrt(2*pi).*exp(-(X.^2/2)-(Y.^2/2)));
surf(X,Y,z);
shading interp
axis tight