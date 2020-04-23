clf;clear all;
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

% Generate transmission signal
rect = zeros(1,N)
rect(mid-half+1:mid+half) = 1
m_1 = conv(rect,rect)

rect = zeros(1,2*N-1)
rect(mid-half+1:mid+half) = 1
m_2 = rect

% Go through the channel
transmission_signal = m_1.*carrier_cos_small+m_2.*carrier_sin_small


% demodulation process
recieve_signal = transmission_signal .* carrier_cos_small
recieve_signal_2 = transmission_signal .* carrier_sin_small

% Go through Low Pass Filter
pole= 1
zero = poly([cos(fc_small*4*pi)+i*sin(fc_small*4*pi),cos(fc_small*4*pi)-i*sin(fc_small*4*pi)])

[h,w] = freqz(zero,pole,'whole',2001)
DC_gain = 10.^(20*log10(abs(h))./20)
recieve_m_1 = filter(zero,pole,recieve_signal)./DC_gain(1)
recieve_m_2 = filter(zero,pole,recieve_signal_2)./DC_gain(1)

subplot(6,1,1);plot(m_1,'.-');title('information 1');grid on;
subplot(6,1,2);plot(m_2,'.-');title('information 2');grid on;
% conduct a modulation
subplot(6,1,3);plot(m_1.*carrier_cos_small,'.-');title('information 1 modulation');grid on;
subplot(6,1,4);plot(abs(fft(transmission_signal)),'.-');title('frequency transmission signal');grid on;
% go through the channel
subplot(6,1,5);plot(recieve_m_1,'.-');title('recieve signal 1');grid on;
subplot(6,1,6);plot(recieve_m_2,'.-');title('recieve signal 2');grid on;
