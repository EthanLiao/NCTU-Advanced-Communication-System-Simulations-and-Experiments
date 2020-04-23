clf;clear all;
fc_small = 1/4
N = 256
mid = ceil(N/2)
half = 128
time = [1:2*N-1]

% modulation process
carrier_cos_small = sqrt(2)*cos(2*pi*fc_small*time)
carrier_sin_small = -sqrt(2)*sin(2*pi*fc_small*time)

% Generate signal
rect = zeros(1,N)
rect(mid-half+1:mid+half) = 1
m_1 = conv(rect,rect)

rect = zeros(1,2*N-1)
rect(mid-half+1:mid+half) = 1
m_2 = rect

% Go through the channel
transmission_signal = m_1 .* carrier_cos_small + m_2 .* carrier_sin_small

% Design a Low Pass Filter
zero= poly([cos(4*fc_small*pi)+i*sin(4*fc_small*pi),cos(4*fc_small*pi)-i*sin(4*fc_small*pi)])
pole= 1

[h,w] = freqz(zero,pole,'whole',2001)
DC_gain = 10.^(20*log10(abs(h))./20)


% conduct a demodulation and filtering
recieve_m_1 = filter(zero,pole,transmission_signal.*carrier_cos_small)./ DC_gain(1)
recieve_m_2 = filter(zero,pole,transmission_signal.*carrier_sin_small)./ DC_gain(1)

% simulate a delay for signal
recieve_m_1_delay = delay(recieve_m_1,30)
recieve_m_2_delay = delay(recieve_m_2,30)

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


function delay_signal = delay(signal,delay)
  delta = zeros(1,delay*2)
  delta(delay) = 1
  delay_signal = conv(delta,signal)
end
