clf;clear all;
carrier_freq = 1/4
N = 50
SNR_db = 20
t =  [1:2*N-1]

% Triangular Wave
half = 10
mid = ceil(N/2)
rect_tri = zeros(1,N)
rect_tri(mid-half:mid+half) = 1
triPulse = conv(rect_tri,rect_tri)

% Rectagular Wave
rec_len = 2*N-1
mid = ceil(rec_len/2)
rect = zeros(1,rec_len)
rect(mid-half:mid+half) = 1

%%%%%%%%%%%%%%%%%%% Real System %%%%%%%%%%%%%%%%%%%%%%%%

carrier_cos = sqrt(2)*cos(2*pi*carrier_freq*t)
carrier_sin = -sqrt(2)*sin(2*pi*carrier_freq*t)

% modulation
orig_sig_1 = rect .* carrier_cos
orig_sig_2 = triPulse .* carrier_sin

tran_sig = orig_sig_1 + orig_sig_2

% AWGN Channel
awgn_sig = add_awgn_noise(tran_sig,SNR_db)

% -------------Reciever------------------
% demodulation
rcv_sig_1 = awgn_sig .* carrier_cos
rcv_sig_2 = awgn_sig .* carrier_sin

% filter it
zero = poly([cos(4*pi*carrier_freq)+i*sin(4*pi*carrier_freq),cos(4*pi*carrier_freq)-i*sin(4*pi*carrier_freq)])
pole = 1

restore_sig_1 = filter(zero,pole,rcv_sig_1)./4
restore_sig_2 = filter(zero,pole,rcv_sig_2)./4

subplot(2,2,1);plot(rect)
subplot(2,2,2);plot(restore_sig_1)
subplot(2,2,3);plot(triPulse)
subplot(2,2,4);plot(restore_sig_2)

%%%%%%%%%%%%%% Equal System %%%%%%%%%%%%%%%

complex_carrier = sqrt(2)*(cos(2*pi*carrier_freq*t) + i*sin(2*pi*carrier_freq*t))
% -----------modulation-------------------
tran_sig = real((rect+i*triPulse) .* complex_carrier)

% AWGN Channel
awgn_sig = add_awgn_noise(tran_sig,SNR_db)

% ----------demodulation-------------------
demod_sig = awgn_sig .* conj(complex_carrier)

zero = poly([cos(4*pi*carrier_freq)+sin(4*pi*carrier_freq),cos(4*pi*carrier_freq)-sin(4*pi*carrier_freq)])
pole = 1
rcv_sig = filter(zero,pole,demod_sig)

restore_sig_1 = real(rcv_sig)
restore_sig_2 = imag(rcv_sig)

figure()
subplot(2,2,1);plot(rect)
subplot(2,2,2);plot(restore_sig_1)
subplot(2,2,3);plot(triPulse)
subplot(2,2,4);plot(restore_sig_2)

function y = add_awgn_noise(x,SNR_db)
  rng('default')
  L = length(x)
  snr = 10^(SNR_db/10)
  SYME = sum(abs(x).^2)/L
  N0 = SYME / snr
  if isreal(x)
    n = sqrt(N0) * randn([1,L])
  else
    n = N0/2 * (randn(1,L)+i*randn(1,L))
  end
  y = x+n
end
