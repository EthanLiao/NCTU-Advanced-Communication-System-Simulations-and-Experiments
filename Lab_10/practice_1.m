clf;clear all;close all;
fc_small = 1/4
N = 256
mid = ceil(N/2)
half = 128
time = [1:2*N-1]
imbalance_phi = 3
g = 3
% modulation process
carrier_cos_small = sqrt(2)*cos(2*pi*fc_small*time)
carrier_sin_small = sqrt(2)*sin(2*pi*fc_small*time)
carrier_sin_small_im = sqrt(2)*sin(2*pi*fc_small*time+imbalance_phi)*g

m_1 = zeros(1, N)
m_1(mid-half+1 : mid+half) = 1

tri_pulse = conv(m_1,m_1)
m_1 = tri_pulse
m_2 = zeros(1, 2*N-1)
m_2(mid-half+1 : mid+half) = 1

% Go through the channel
complex_carrier = carrier_cos_small + i*carrier_sin_small_im
transmission_signal = real((m_1 + i*m_2) .* complex_carrier)


% conduct a demodulation
recieve_m = transmission_signal.*conj(carrier_cos_small + i*carrier_sin_small)

% Design a Low Pass Filter
pole= 1
zero= poly([cos(4*fc_small*pi)+i*sin(4*fc_small*pi),cos(4*fc_small*pi)-i*sin(4*fc_small*pi)])

[h,w] = freqz(zero,pole,'whole',2001)
DC_gain = 10.^(20*log10(abs(h))./20)
recieve_m_tmp = filter(zero,pole,recieve_m)./ DC_gain(1)
recieve_m_1 = real(recieve_m_tmp)
recieve_m_2 = imag(recieve_m_tmp)

subplot(3,2,1);plot(m_1,'.-');title('information 1');grid on;
subplot(3,2,2);plot(m_2,'.-');title('information 2');grid on;
% conduct a modulation
subplot(3,2,3);plot(abs(fft(transmission_signal)),'.-');title('frequency transmission signal');grid on;
% go through the channel
subplot(3,2,4);plot(abs(fft(recieve_m)),'.-');title('frequency transmission signal');grid on;
subplot(3,2,5);plot(recieve_m_1,'.-');title('recieve information 1');grid on;
subplot(3,2,6);plot(recieve_m_2,'.-');title('recieve information 2');grid on;
