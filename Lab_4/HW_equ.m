clf
fc_small = 1/4
N = 256
mid = ceil(N/2)
half = 128
time = [1:2*N-1]

% modulation process
carrier_cos_small = sqrt(2)*cos(2*pi*fc_small*time)
carrier_sin_small = sqrt(2)*sin(2*pi*fc_small*time)

m_1 = zeros(1, N)
m_1(mid-half+1 : mid+half) = 1



tri_pulse = conv(m_1,m_1)
m_1 = tri_pulse
m_2 = zeros(1, 2*N-1)
m_2(mid-half+1 : mid+half) = 1

% Go through the channel
complex_carrier = carrier_cos_small + i*carrier_sin_small
transmission_signal = real((m_1 + i*m_2) .* complex_carrier)




% conduct a demodulation
recieve_m = transmission_signal.*conj(complex_carrier)

% Design a Low Pass Filter
pole= 1
zero= poly([cos(4*fc_small*pi)+i*sin(4*fc_small*pi),cos(4*fc_small*pi)-i*sin(4*fc_small*pi)])

[h,w] = freqz(zero,pole,'whole',2001)
DC_gain = 10.^(20*log10(abs(h))./20)
recieve_m_tmp = filter(zero,pole,recieve_m)./ DC_gain(1)
recieve_m_1 = real(recieve_m_tmp)
recieve_m_2 = imag(recieve_m_tmp)

% Design a delay in channel
delay = 6
transmission_signal= ifft(fft(transmission_signal)*exp(i*2*pi*fc_small*delay))
% conduct a demodulation
recieve_m = transmission_signal.*conj(complex_carrier)
recieve_m_tmp = filter(zero,pole,recieve_m)./ DC_gain(1)
recieve_m_1_delay = real(recieve_m_tmp)
recieve_m_2_delay = imag(recieve_m_tmp)

subplot(5,2,1);plot(m_1,'.-');title('information 1');grid on;
subplot(5,2,2);plot(m_2,'.-');title('information 2');grid on;
% conduct a modulation
subplot(5,2,3);plot(m_1.*carrier_cos_small,'.-');title('information 1 modulation');grid on;
subplot(5,2,4);plot(m_2.*carrier_cos_small,'.-');title('information 2 modulation');grid on;
% go through the channel
subplot(5,2,5);plot(recieve_m_1,'.-');title('recieve information 1');grid on;
subplot(5,2,6);plot(recieve_m_2,'.-');title('recieve information 2');grid on;
subplot(5,2,7);plot(recieve_m_1_delay,'.-');title('recieve delay information 1');grid on;
subplot(5,2,8);plot(recieve_m_2_delay,'.-');title('recieve delay information 2');grid on;
subplot(5,2,9);plot(w/pi,20*log10(abs(h)),'.-');title('Frequency domain low pass filter');grid on;
