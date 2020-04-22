Ts   = 1e-6; % Symbol time (sec)
span = 6;    % Filter span in symbols

a = Ts*[.5, .75, 1, 2];
B = sqrt(log(2)/2)./(a);

t = linspace(-span*Ts/2,span*Ts/2,1000)';
hg = zeros(length(t),length(a));
for k = 1:length(a)
  hg(:,k) = sqrt(pi)/a(k)*exp(-(pi*t/a(k)).^2);
end

figure(1)
plot(t/Ts,hg)
title({'Impulse response of a continuous-time Gaussian filter';...
  'for various bandwidths'});
xlabel('Normalized time (t/Ts)')
ylabel('Amplitude')
legend(sprintf('a = %g*Ts',a(1)/Ts),sprintf('a = %g*Ts',a(2)/Ts),...
  sprintf('a = %g*Ts',a(3)/Ts),sprintf('a = %g*Ts',a(4)/Ts))
grid on;

% Frequency Response for Continuous-time Gaussian Filter

f = linspace(0,32e6,10000)';
Hideal = zeros(length(f),length(a));
for k = 1:length(a)
  Hideal(:,k) = exp(-a(k)^2*f.^2);
end

figure(2)
plot(f,20*log10(Hideal))
titleStr = {'Ideal magnitude response for a continuous-time ';...
  'Gaussian filter for various bandwidths'};
title(titleStr);
legend(sprintf('B = %g',B(1)),sprintf('B = %g',B(2)),...
  sprintf('B = %g',B(3)),sprintf('B = %g',B(4)))
hold on
for k = 1:length(a)
  plot(B,20*log10(exp(-a.^2.*B.^2)),'ro','HandleVisibility','off')
end

axis([0 5*max(B) -50 5])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on;
