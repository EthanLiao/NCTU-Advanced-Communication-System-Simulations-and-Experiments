clf;clear all;
N = 3
fc = 1/6
t = 0:0.001:N / fc
s0 = sin(2*pi*fc*t)
s1 = sin(2*pi*fc*t+pi)
% signal = [s0,s1]
signal = [1,1,1,-1,1,-1,1,1,1]
t_vec = [t t(end)+t]
load('HW_filter')
load('TEST')
symbol_rate = 1*10^6
DAC_sampling_factor = (4*10^6)/symbol_rate
DMA_sampling_factor = (32*10^6)/(4*10^6)
carrier_frequency = 8*10^6
zero_freq = 1/40
pole_freq = 1

srrc = srrc_pulse(4, 11, 5, 0.2);

zero = [0.028  0.053 0.071  0.053 0.028]
pole = [1.000 -2.026 2.148 -1.159 0.279]



% modulatoin part
t_Digital_filtered_signal = conv(DAC(signal,DAC_sampling_factor),srrc)
t_DMA_filtered_signal = filter(TEST,DAC(t_Digital_filtered_signal,DMA_sampling_factor))

% modulatoin with carrier
t = [0:length(t_DMA_filtered_signal)-1]
carrier = cos(2*pi*carrier_frequency*t)+i*sin(2*pi*carrier_frequency*t)
signal_with_carrier = real(t_DMA_filtered_signal .* carrier)

% demodulation
demodulation_signal = signal_with_carrier .* conj(carrier)
r_DMA_filtered_signal = ADC(filter(TEST,demodulation_signal),DMA_sampling_factor)
r_Digital_filtered_signal = ADC(conv(r_DMA_filtered_signal,srrc),DAC_sampling_factor)

subplot(4,2,1);stem(signal);title('origin signal')
subplot(4,2,2);plot(abs(fft(signal)));title('origin signal f domain')
subplot(4,2,3);stem(DAC(t_Digital_filtered_signal,DMA_sampling_factor));title('transmit signal')
subplot(4,2,4);plot(abs(fft(DAC(t_Digital_filtered_signal,DMA_sampling_factor))));title('transmit signal f domain')
subplot(4,2,5);stem(real(r_Digital_filtered_signal)./max(abs(real(r_Digital_filtered_signal))));title('receive signal')
subplot(4,2,6);plot(abs(fft(real(r_Digital_filtered_signal))));title('receive signal f domain')
subplot(4,2,7);stem(real(r_Digital_filtered_signal)./max(abs(real(r_Digital_filtered_signal))));hold on;stem(shift(signal,3));title('comparison two signal')


function shift_signal = shift(signal,shift)
  shift_signal = [zeros(1,shift-1) signal]
end

% subplot(3,2,5);stem(r_Digital_filtered_signal./max(abs(r_Digital_filtered_signal)));title('receive signal')
% subplot(3,2,6);plot(abs(fft(r_Digital_filtered_signal./max(r_Digital_filtered_signal))));title('receive signal f domain')
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
