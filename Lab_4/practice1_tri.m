fc_small = 1/64
fc_large = 1/4
N = 256
mid = ceil(N/2)
half = 128
time = [1:2*N-1]

% modulation process
carrier_cos_small = sqrt(2)*cos(2*pi*fc_small*time)
carrier_cos_large = sqrt(2)*cos(2*pi*fc_large*time)
carrier_sin_small = -sqrt(2)*sin(2*pi*fc_small*time)
carrier_sin_large = -sqrt(2)*sin(2*pi*fc_large*time)

m_1 = zeros(1, N)
m_1(mid-half+1 : mid+half) = 1
m_2 = zeros(1, N)
m_2(mid-half+1 : mid+half) = 1

tri_pulse = conv(m_1,m_1)
m_1 = tri_pulse
m_2 = zeros(1, 2*N-1)
m_2(mid-half+1 : mid+half) = 1

% Go through the channel
transmission_signal = m_1.*carrier_cos_large+m_2.*carrier_sin_large


% demodulation process
recieve_signal = transmission_signal .* carrier_cos_large
recieve_signal_2 = transmission_signal .* carrier_sin_large
% Go through Low Pass Filter
pole= 1
zero = poly([cos(fc_large*4*pi)+i*sin(fc_large*4*pi),cos(fc_large*4*pi)-i*sin(fc_large*4*pi)])

[h,w] = freqz(zero,pole,'whole',2001)
DC_gain = 10.^(20*log10(abs(h))./20)
recieve_m_1 = filter(zero,pole,recieve_signal)./DC_gain(1)
recieve_m_2 = filter(zero,pole,recieve_signal_2)./DC_gain(1)

subplot(6,1,1);plot(m_1,'.-');title('information');grid on;
subplot(6,1,2);plot(m_2,'.-');title('cosine carrier');grid on;
% conduct a modulation
subplot(6,1,3);plot(m_1.*carrier_cos_small,'.-');title('modulation');grid on;
subplot(6,1,4);plot(abs(fft(transmission_signal)),'.-');title('frequency transmission signal');grid on;
% go through the channel
subplot(6,1,5);plot(recieve_m_2,'.-');title('recieve signal 2');grid on;
% subplot(10,1,6);plot(recieve_signal,'.-');title('recieve signal');grid on;
subplot(6,1,6);plot(recieve_m_1,'.-');title('recieve signal 1');grid on;
