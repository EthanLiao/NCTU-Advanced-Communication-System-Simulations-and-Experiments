% plot rc pulse
% rc = rcosdesign(0.55,5,4)
rc = rcro_pulse(50,10,2,0.01)
subplot(3,2,1);plot(rc);title('raised cosine');grid on;
subplot(3,2,2);plot(abs(fft(rc)));title('raised cosine f domain');grid on;
% fvtool(rc,'Analysis','impulse');;grid on;
% plot sqrc pulse
Fs = 10; %Sampling frequency
T = 1/Fs; %Symbol period
L = 1000; %Number of points
t = (0:L-1)*T; %Time vector
phi1 = srrc_pulse(1, 1/10, 4, 0); %SRRC pulse
NFFT = 2^nextpow2(L);
Y = fft(phi1, NFFT)./L;
Y2 = fftshift(fft(phi1, NFFT)./L);
f = Fs/2*linspace(0,1,NFFT/2+1);
f2 = ((0:NFFT-1) -ceil((NFFT-1)/2))/NFFT/T;
subplot(3,2,3);plot(phi1);title('square root raise cosine');grid on;
subplot(3,2,4);plot(f2, 2*abs(Y2));title('square root raise cosine f domain');grid on;
subplot(3,2,5);plot(conv(phi1,phi1));title('convolution of two srrc');grid on;
subplot(3,2,6);plot(abs(fft(conv(phi1,phi1))));title('convolution of two srrc f domain');grid on;


function [ pulse_vector, time_vector ] = rcro_pulse( k_t, bit_period, samples_per_bit, r )
% -----------
% k_t = 50
% bit_period = T = 10
% samples_per_bit = sampling_factor = 1
% r = a = 0.5
%rcro_pulse Generate Raised Cosine Rolloff pulse vector
    f0 = 1/(bit_period*2);
    f_delta = r*f0;
    stop_time = -k_t*bit_period;
    start_time= k_t*bit_period;
    time_vector = linspace(stop_time, start_time, samples_per_bit*2*k_t);

    pulse_vector_1 = ((2*f0*sin(2*pi*f0.*time_vector)  ./ ...
                        (2*pi*f0.*time_vector)));
    % Separate the second half of the equation to handle small roundoff
    % errors and 0/0
    pulse_vector_2 = cos(2*pi*f_delta*time_vector);
    pulse_vector_3 = (1-(4*f_delta.*time_vector).^2);
    pulse_vector_2(abs(pulse_vector_2) < 5e-14) = 0;
    pulse_vector_3(abs(pulse_vector_3) < 5e-14) = 0;
    pulse_vector_2 = pulse_vector_2 ./ pulse_vector_3;
    % Apply L'Hopital's rule
    pulse_vector_2(find(isnan(pulse_vector_2))) = ...
                        sin(2*pi*f_delta*time_vector(find(isnan(pulse_vector_2)))) ...
                        * pi ./ ...
                        (16*f_delta*time_vector(find(isnan(pulse_vector_2))));
    % Replace potential NaN with one, mainly at t=0.
    pulse_vector = pulse_vector_1 .* pulse_vector_2;
    pulse_vector(find(isnan(pulse_vector))) = 2*f0;
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
