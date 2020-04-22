clf;clear all
N = 3
fc = 6
t = 0:0.001:N / (2*pi)
s0 = sin(2*pi*fc*t)
s1 = sin(2*pi*fc*t+pi)
signal = s0

% delayed signal
delay_factor_fast = 10 ; delay_factor_slow = 20
channel_gain_fast = 0.8; channel_gain_slow = 0.9
width = delay_factor_fast*2-1
delta = zeros(1,width)
delta(delay_factor_fast) = 1
delayed_signal_fast = conv(s0,delta)
width = delay_factor_slow*2-1
delta = zeros(1,width)
delta(delay_factor_slow) = 1
delayed_signal_slow = conv(s0,delta)



fast_signal = channel_gain_fast.*signal_go_through_channel(delayed_signal_fast)
slow_signal = channel_gain_slow.*signal_go_through_channel(delayed_signal_fast)
ISI_signal = [fast_signal;zeros(1,length(slow_signal)-length(fast_signal))]+slow_signal

% Use equalzer to equalize signal
zero = [channel_gain_fast , channel_gain_slow]
pole = 1
equalized_signal = filter(zero,pole,ISI_signal)

subplot(3,1,1);plot(signal);title('Transmitted Signal');grid on;
subplot(3,1,2);plot(ISI_signal);title('Without Equalize Recieved Signal(ISI Signal)');grid on;
subplot(3,1,3);plot(equalized_signal);title('Equalized Recieved Signal');grid on;


function r_Digital_filtered_signal = signal_go_through_channel(signal)
  symbol_rate = 1*10^6
  DAC_sampling_factor = 16*10^6/symbol_rate
  DMA_sampling_factor = 32*10^6/symbol_rate
  carrier_frequency = 8*10^6
  srrc = srrc_pulse(1/symbol_rate, 1/10, 4, 0);
  AWGN_SNR = 5

  % modulatoin part
  t_Digital_filtered_signal = conv(DAC(signal,DAC_sampling_factor),srrc)
  t_DMA_filtered_signal = conv(DAC(t_Digital_filtered_signal,DMA_sampling_factor),srrc)

  % modulatoin with carrier
  t = [0:length(t_DMA_filtered_signal)-1]
  carrier = cos(2*pi*carrier_frequency*t)+i*sin(2*pi*carrier_frequency*t)
  signal_with_carrier = real(t_DMA_filtered_signal .* carrier)
  signal_with_carrier = add_awgn_noise(signal_with_carrier,AWGN_SNR)
  % demodulation
  demodulation_signal = signal_with_carrier .* conj(carrier)
  r_DMA_filtered_signal = ADC(conv(demodulation_signal,srrc),DMA_sampling_factor)
  r_Digital_filtered_signal = real(ADC(conv(r_DMA_filtered_signal,srrc),DAC_sampling_factor))
end

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

function y = add_awgn_noise(x,SNR_dB)
 %y=awgn_noise(x,SNR) adds AWGN noise vector to signal 'x' to generate a
 %resulting signal vector y of specified SNR in dB
 rng('default');%set the random generator seed to default (for comparison only)
 L=length(x);
 SNR = 10^(SNR_dB/10); %SNR to linear scale
 Esym=sum(abs(x).^2)/(L); %Calculate actual symbol energy
 N0=Esym/SNR; %Find the noise spectral density
 if(isreal(x)),
 noiseSigma = sqrt(N0);%Standard deviation for AWGN Noise when x is real
 n = noiseSigma*randn(1,L);%computed noise
 else
 noiseSigma=sqrt(N0/2);%Standard deviation for AWGN Noise when x is complex
 n = noiseSigma*(randn(1,L)+1i*randn(1,L));%computed noise
 end
 y = x + n; %received signal
end
