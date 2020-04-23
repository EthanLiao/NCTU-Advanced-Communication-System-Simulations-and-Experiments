clf;clear all;
% plot bpsk signal
load('./filter/Pulse_Shapping_Filter')
load('./filter/FIR_distortion')

signal = [-1,1,-1,1]
delayed_signal = delay(signal,0)

% sout = filter(FIR_distortion,signal)
transmit_signal = filter(Pulse_Shapping_Filter,DAC(signal))
recieve_signal = ADC(filter(Pulse_Shapping_Filter,transmit_signal))

stem(delayed_signal)
hold on;stem(recieve_signal./0.056);grid on;title('Pratical IIR Pulse Shapping signal')

function delay_sig = delay(sig,damt)
  delay_sig = [zeros(1,damt) sig]
end

function ADC_sig = ADC(origin_signal)
down = 32
ADC_sig = origin_signal([1:down:length(origin_signal)])
end

function DAC_sig = DAC(origin_signal)
up = 32
DAC_sig = zeros(1,up*length(origin_signal))
DAC_sig([1:up:length(DAC_sig)]) = origin_signal
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
