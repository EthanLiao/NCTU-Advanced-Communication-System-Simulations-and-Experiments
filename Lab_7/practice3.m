clf;clear all;
% plot bpsk signal
load('Pulse_Shapping_Filter')
load('FIR_distortion')
delta_25 = zeros(1,100)
delta_25(25) = 1
[signal,bit] = bpskd([1 0],2)
delayed_signal = conv(signal,delta_25)
[h,w] = freqz(FIR_distortion,'whole',2001)

FIR_DC_gain = 10.^(20*log10(abs(h))./20)
sout = filter(FIR_distortion,signal)./FIR_DC_gain(1)

[h,w] = freqz(Pulse_Shapping_Filter,'whole',2001)
Pulse_Shapping_gain = 10.^(20*log10(abs(h))./20)


transmit_signal = filter(Pulse_Shapping_Filter,DAC(sout))./Pulse_Shapping_gain(1)
recieve_signal = ADC(filter(Pulse_Shapping_Filter,transmit_signal)./Pulse_Shapping_gain(1))./0.000838

plot(delayed_signal)
hold on;plot(recieve_signal);grid on;title('Pratical IIR Pulse Shapping signal')

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


function [bpsk,bit]=bpskd(g,f)
if nargin > 2
    error('Too many input arguments');
elseif nargin==1
    f=1;
end
if f<1;
    error('Frequency must be bigger than 1');
end
t=0:2*pi/99:2*pi;
cp=[];sp=[];
mod=[];mod1=[];bit=[];
for n=1:length(g);
    if g(n)==0;
        die=-ones(1,100);   %Modulante
        se=zeros(1,100);    %Señal
    else g(n)==1;
        die=ones(1,100);    %Modulante
        se=ones(1,100);     %Señal
    end
    c=sin(f*t);
    cp=[cp die];
    mod=[mod c];
    bit=[bit se];
end
bpsk=cp.*mod;
end
