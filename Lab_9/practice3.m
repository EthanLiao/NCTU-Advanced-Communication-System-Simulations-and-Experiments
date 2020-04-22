N = 3
fc = 1/6
t = 0:0.001:N / fc
s0 = sin(2*pi*fc*t)
s1 = sin(2*pi*fc*t+pi)
signal = [s0,s1]
t_vec = [t t(end)+t]

symbol_rate = 1*10^6
IF_frequency = 4*10^6
ADC_sampling_factor = 16*10^6/symbol_rate
DMA_sampling_factor = 32*10^6/symbol_rate
carrier_frequency = 8*10^6
srrc = srrc_pulse(1/symbol_rate, 1/10, 4, 0);
srrc_dc_gain = abs(fft(srrc))

% modulatoin part
carrier = cos(2*pi*(carrier_frequency-IF_frequency)*t_vec)
signal_with_carrier = signal .* carrier

t_DMA_filtered_signal = DAC(conv(signal_with_carrier,srrc)./srrc_dc_gain(1),DMA_sampling_factor)
% t_DMA_filtered_signal = conv(ADC(t_Digital_filtered_signal,DMA_sampling_factor),srrc)

% de modulatoin with carrier
t = [0:length(t_DMA_filtered_signal)-1]
carrier = cos(2*pi*IF_frequency*t)-i*sin(2*pi*IF_frequency*t)
demodulation_signal = t_DMA_filtered_signal .* carrier

% demodulation
r_Digital_filtered_signal = ADC(conv(demodulation_signal,srrc)./srrc_dc_gain(1),ADC_sampling_factor)

subplot(2,1,1);plot(signal);title('BPSK Signal');grid on;
subplot(2,1,2);plot(real(r_Digital_filtered_signal));title('Recieved Signal');grid on;

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
