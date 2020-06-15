clf;clear all;close all;

N = 1000;
sig_1 = randi([0,1],1,N);
sig_1((sig_1==0)) = -1;

sig_2 = randi([0,1],1,N);
sig_2((sig_2==0)) = -1;

sig = [sig_1 ; sig_2];


% 2x2 System
H = randn(2,2);
H = sqrt(1/2) * H +1j * sqrt(1/2) * H;
mimo_rcv = H*sig;


SNR_DB = 5;
% recieve signal
mimo_rcv_1 = add_awgn_noise(mimo_rcv(1,:), SNR_DB);
mimo_rcv_2 = add_awgn_noise(mimo_rcv(2,:), SNR_DB);


mimo_rcv = [mimo_rcv_1 ; mimo_rcv_2];

% ------ Zero Forcing -------------
inv_H = inv(H);
rcv_beam = inv_H * mimo_rcv;

% signal detect
rcv_beam = real(rcv_beam);
rcv_beam(rcv_beam>0) = 1;
rcv_beam(rcv_beam<0) = -1;

% symbol error rate
zf_symbol_err = 1-mean(sig==rcv_beam, 'all')
zf_snr = SNR(sig,rcv_beam);

% Zero forcing
subplot(2,2,1);stem(sig_1);title("ZF transmission signal branch 1");
subplot(2,2,2);stem(sig_2);title("ZF transmission signal branch 2");
subplot(2,2,3);stem(rcv_beam(1,:));title("ZF recieved signal branch 1");
subplot(2,2,4);stem(rcv_beam(2,:));title("ZF recieved signal branch 2");

% ------------------MMSE Process-----------------

% MMSE Dectector
% calculate Signal to noise power ratio : sigm

sigm = 10^(SNR_DB/10);

W = inv(H'*H + 1/sigm * eye(2))*H';
mmse_rcv_beam = W * mimo_rcv;

% signal detect
mmse_rcv_beam = real(mmse_rcv_beam);
% rcv_beam = rcv_beam / max(rcv_beam,[],'all');
mmse_rcv_beam(mmse_rcv_beam>0) = 1;
mmse_rcv_beam(mmse_rcv_beam<0) = -1;
%
% % symbol error rate
mmse_symbol_err = 1-mean(sig==mmse_rcv_beam, 'all')
mmse_snr = SNR(sig, mmse_rcv_beam);

% MMSE
figure();
subplot(2,2,1);stem(sig_1);title("MMSE transmission signal branch 1");
subplot(2,2,2);stem(sig_2);title("MMSE transmission signal branch 2");
subplot(2,2,3);stem(mmse_rcv_beam(1,:));title("MMSE recieved signal branch 1");
subplot(2,2,4);stem(mmse_rcv_beam(2,:));title("MMSE recieved signal branch 2");

function y = add_awgn_noise(x,SNR_DB)
  L = length(x);
  % calculate symbol energy
  SNR_N = 10^(SNR_DB/10); % SNR enery to linear scale
  SYME = sum(abs(x).^2) / L;
  N0 = SYME / SNR_N;      % Noise spectral Density
  if isreal(x)
    n = sqrt(N0) * randn(1,L);
  else
    n = sqrt(N0/2) * (randn(1,L)+i*randn(1,L));
  end
  y = x + n;
end

function snr = SNR(x,y)
  snr = 10*log10(mean(x.^2, 'all') / mean((y-x).^2, 'all'));
end
