clf;clear all;
N = 3
fc = 1/6
t = (0:0.001:N)/ fc
symbol_rate = 0.1
% s0 = sin(2*pi*fc*t)
% s1 = sin(2*pi*fc*t+pi)
% signal = [s0,s1]
% t_vec = [t t(end)+t]
signal = [-1,1,-1,1,1,1]
load('./filter/practice_1_filter')
srrc = srrc_pulse(4 , 4, 0.2);
Digital_filtered_signal = conv(DAC(signal,4),srrc)
DMA_filtered_signal = DAC(Digital_filtered_signal,16)
DMA_analog_signal = filter(practice_1_filter,DMA_filtered_signal)
subplot(2,1,1);stem(signal);title('BPSK Signal');grid on;
subplot(2,1,2);stem(DMA_analog_signal);title('Practical DAC');grid on;


function DAC_sig = DAC(origin_signal,up_factor)
DAC_sig = zeros(1,up_factor*length(origin_signal))
DAC_sig([1:up_factor:length(DAC_sig)]) = origin_signal
end


function [phi, t] = srrc_pulse(T, A, a)
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
t = [-A*T:A*T] + 10^(-8); % in order to avoid division by zero problems at t=0.
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
