clf;clear all; close all;
N = 3
fc = 1/6
t = 0:0.001:N / fc
signal = [-1,1,1,-1,1]
load('IIR_filter')

symbol_rate = 1*10^6
carrier_frequency = 8*10^6
srrc = srrc_pulse(32, 11, 5, 0.3);

% modulatoin part
t_Digital_filtered_signal = conv(DAC(signal,16),srrc,'same')
t_DMA_filtered_signal = filter(IIR_filter,DAC(t_Digital_filtered_signal,2))

% modulatoin with carrier
t = [0:length(t_DMA_filtered_signal)-1]
carrier = cos(2*pi*carrier_frequency*t)+i*sin(2*pi*carrier_frequency*t)
signal_with_carrier = real(t_DMA_filtered_signal .* carrier)

% demodulation
demodulation_signal = signal_with_carrier .* conj(carrier)
r_Digital_filtered_signal = real(ADC(filter(IIR_filter,demodulation_signal),32))

recieve_signal = r_Digital_filtered_signal./max(abs(r_Digital_filtered_signal))
r_signal = (recieve_signal>0)+(-(recieve_signal<0))

subplot(2,1,1);stem(delay(signal,0));title('BPSK Signal');grid on;
subplot(2,1,2);stem(r_signal);title('Recieved Signal');grid on;

function delay_sig = delay(orig_sig,damt)
  delay_sig = [zeros(1,damt) orig_sig]
end
function DAC_sig = DAC(origin_signal,up_factor)
DAC_sig = zeros(1,up_factor*length(origin_signal))
DAC_sig([1:up_factor:end]) = origin_signal
end

function ADC_sig = ADC(origin_signal,down_factor)
ADC_sig = zeros(1,length(origin_signal))
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
