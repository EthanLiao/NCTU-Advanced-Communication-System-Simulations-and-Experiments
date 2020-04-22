clf;clear all;
N = 3
fc = 1/6
NOISE_DB = 100
t = 0:0.001:N / fc
s0 = sin(2*pi*fc*t)
s1 = sin(2*pi*fc*t+pi)
% signal = [s0,s1]
signal = [1,-1,1,-1,-1,1,1,1,1]
t_vec = [t t(end)+t]
srrc_4 = srrc_pulse(4,5, 0.2)
srrc_16 = srrc_pulse(16,5, 0.2)
delta = zeros(1,40)
delta(12) = 1
load('practice_2_filter')
t_Digital_filtered_signal = conv(DAC(signal,4),srrc_4)
t_DMA_filtered_signal = filter(practice_2_filter,DAC(t_Digital_filtered_signal,4))
channel_signal = add_awgn_noise(t_DMA_filtered_signal,NOISE_DB)
r_DMA_filtered_signal = ADC(filter(practice_2_filter,channel_signal),4)
r_Digital_filtered_signal = ADC(conv(r_DMA_filtered_signal,srrc_4),4)

delay_signal = conv(signal,delta)
stem(delay_signal./max(delay_signal));title('BPSK Signal');grid on;
% subplot(3,1,2);plot(t_DMA_filtered_signal);title('Practical DAC');grid on;
hold on;stem(r_Digital_filtered_signal./max(r_Digital_filtered_signal));grid on;

function DAC_sig = DAC(origin_signal,up_factor)
DAC_sig = zeros(1,up_factor*length(origin_signal))
DAC_sig([1:up_factor:length(DAC_sig)]) = origin_signal
end

function ADC_sig = ADC(origin_signal,down_factor)
ADC_sig = origin_signal([1:down_factor:length(origin_signal)])
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
