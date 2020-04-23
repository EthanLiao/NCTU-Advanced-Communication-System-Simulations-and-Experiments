N = 3
fc = 1/6
t = (0:0.001:N )/ fc
% ------ BPSK Signal-----------
s0 = sin(2*pi*fc*t)
s1 = sin(2*pi*fc*t+pi)
signal = [s0,s1]
t_vec = [t t(end)+t]
% ------------------------------
% ------- QPSK Signal-----------
data=[0 1 0 1]; % information
data_NZR=2*data-1; % Data Represented at NZR form for QPSK modulation
s_p_data=reshape(data_NZR,2,length(data)/2);  % S/P convertion of data
br=10.^6; %Let us transmission bit rate  1000000
f=br; % minimum carrier frequency
T=1/br; % bit duration
t=T/99:T/99:T; % Time vector for one bit information
y=[];
y_in=[];
y_qd=[];
for(i=1:length(data)/2)
    y1=s_p_data(1,i)*cos(2*pi*f*t); % inphase component
    y2=s_p_data(2,i)*sin(2*pi*f*t) ;% Quadrature component
    y_in=[y_in y1]; %inphase signal vector
    y_qd=[y_qd y2]; %quadrature signal vector
    y=[y y1+y2]; % modulated signal vector
end
Tx_sig=y; % transmitting signal after modulation
tt=T/99:T/99:(T*length(data))/2;
% -----------------

symbol_rate = 1*10^6
DAC_sampling_factor = 64
DMA_sampling_factor = 128
carrier_frequency = 8*10^6
srrc = srrc_pulse(1/symbol_rate, 1/10, 4, 0);

% modulatoin part
t_Digital_filtered_signal = conv(DAC(Tx_sig,DAC_sampling_factor),srrc)
t_DMA_filtered_signal = conv(DAC(t_Digital_filtered_signal,DMA_sampling_factor),srrc)

% modulatoin with carrier
t = [0:length(t_DMA_filtered_signal)-1]
carrier = cos(2*pi*carrier_frequency*t)+i*sin(2*pi*carrier_frequency*t)
signal_with_carrier = real(t_DMA_filtered_signal .* carrier)

% demodulation
demodulation_signal = signal_with_carrier .* conj(carrier)
r_DMA_filtered_signal = ADC(conv(demodulation_signal,srrc),DMA_sampling_factor)
r_Digital_filtered_signal = ADC(conv(r_DMA_filtered_signal,srrc),DAC_sampling_factor)./55.47.*1.4142

subplot(3,2,1);plot(Tx_sig);title('QPSK Signal');grid on;
subplot(3,2,2);plot(abs(fft(Tx_sig)));title('QPSK Signal');grid on;
subplot(3,2,3);plot(t_Digital_filtered_signal);title('QPSK Pulse Shaping Signal');grid on;
subplot(3,2,4);plot(abs(fft(t_Digital_filtered_signal)));title('QPSK Pulse Shaping Signal f domain');grid on;
subplot(3,2,5);plot(real(r_Digital_filtered_signal));title('Recieved Signal');grid on;
subplot(3,2,6);plot(abs(fft(real(r_Digital_filtered_signal))));title('Recieved Signal');grid on;

function DAC_sig = DAC(origin_signal,up_factor)
DAC_sig = zeros(1,up_factor*length(origin_signal))
DAC_sig([1:up_factor:length(DAC_sig)]) = origin_signal
end

function ADC_sig = ADC(origin_signal,down_factor)
ADC_sig = origin_signal([1:down_factor:length(origin_signal)])
end

function [phi, t] = srrc_pulse(T, Ts, A, a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi = srrc_pulse(T, Ts, A, a)                                                 %
% OUTPUT                                                                        %
%      phi: truncated SRRC pulse, with parameter T,                             %
%                 roll-off factor a, and duration 2*A*T                         %
%      t:   time axis of the truncated pulse                                    %
% INPUT                                                                         %
%      T:  Nyquist parameter or symbol period  (real number)                    %
%      Ts: sampling period  (Ts=T/over)                                         %
%                where over is a positive INTEGER called oversampling factor    %
%      A:  half duration of the pulse in symbol periods (positive INTEGER)      %
%      a:  roll-off factor (real number between 0 and 1)                        %
%                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = [-A*T:Ts:A*T] + 10^(-8); % in order to avoid division by zero problems at t=0.
if (a>0 && a<=1)
   num = cos((1+a)*pi*t/T) + sin((1-a)*pi*t/T) ./ (4*a*t/T);
   denom = 1-(4*a*t./T).^2;
   phi = 4*a/(pi*sqrt(T)) * num ./ denom;
elseif (a==0)
   phi = 1/(sqrt(T)) * sin(pi*t/T)./(pi*t/T);
else
    phi = zeros(length(t),1);
    disp('Illegal value of roll-off factor')
    return
end
end
