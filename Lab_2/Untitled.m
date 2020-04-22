clc;
clearvars;
close all;
imtool close all;  % Close all imtool figures.
workspace;
format longg;
format compact;
fontSize = 16;
numberOfElements = 1001;
middleElement = ceil(numberOfElements/2);
x = linspace(-10,10, numberOfElements);
halfWindowWidth = 150;
% Make rectangular pulse.
rect = zeros(1, numberOfElements);
rect(middleElement-halfWindowWidth: middleElement+halfWindowWidth) = 1;

plot(x, rect,'b-', 'LineWidth', 2);
ylim([0 2]);
grid 'on';
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

% Convolve
y = conv(rect, rect,'full');
predictedLength = length(rect)+length(rect)-1
newXAxis = linspace((-10 - halfWindowWidth), (+10 + halfWindowWidth), length(y)) ;

% plot(newXAxis, y, 'b-', 'LineWidth', 2);
grid 'on';

var_x = var(y)
% var_noise = var_x ^ 2 / (10^(SNR/10))
s = 10000;
noise =  10*randn(1,2001)
rect_with_noise = noise+y
fft_y =  fft(rect_with_noise)

subplot(3, 1, 1);
plot(rect_with_noise)
% plot(ifft(win(fft_y)))
subplot(3, 1, 2);


plot(fft_y)

for idx = 1:2001
  if abs(fft_y(idx)) < 2000
    fft_y(idx) = 0
  end
end
% plot(ifft(fft_noise))
subplot(3, 1, 3);
plot(ifft(fft_y))
