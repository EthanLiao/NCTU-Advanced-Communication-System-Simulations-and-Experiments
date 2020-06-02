clf;clear all;close all;

N = 10;
sig_1 = randi([0,1],1,N);
sig_1((sig_1==0)) = -1;
sig_2 = randi([0,1],1,N);
sig_2((sig_2==0)) = -1;
sig = [sig_1 ; sig_2];

SNR_l = [-20:1:20];
zf_snr = zeros(1,length(SNR_l));
mmse_snr = zeros(1,length(SNR_l));

for i=1:length(SNR_l)
  [zf_snr(i), mmse_snr(i)] = MIMO_eval(SNR_l(i), sig);
end


plot(SNR_l,zf_snr);
hold on;
plot(SNR_l,mmse_snr);
legend("Zero forcing SNR", "MMSE SNR")

function  [zf_snr, mmse_snr]= MIMO_eval(SNR_DB, sig)
  % 2x2 System
  H = randn(2,2);
  H = sqrt(1/2) * H +1j * sqrt(1/2) * H;
  mimo_rcv = H*sig;
  % recieve signal
  mimo_rcv_1 = add_awgn_noise(mimo_rcv(1,:), SNR_DB);
  mimo_rcv_2 = add_awgn_noise(mimo_rcv(2,:), SNR_DB);

  mimo_rcv = [mimo_rcv_1 ; mimo_rcv_2];

  % ------ Zero Forcing -------------
  inv_H = inv(H);
  rcv_beam = inv_H * mimo_rcv;
  % signal detect
  rcv_beam = real(rcv_beam);
  % rcv_beam = rcv_beam / max(rcv_beam,[],'all');
  % rcv_beam(rcv_beam>0.5) = 1;
  % rcv_beam(rcv_beam<-0.5) = -1;
  % symbol error rate
  % zf_symbol_err = 1-mean(sig==rcv_beam, 'all')
  zf_snr = SNR(sig,rcv_beam);
  % ------------------MMSE Process-----------------

  % MMSE Dectector
  % calculate Signal to noise power ratio : sigm
  sigm = 10^(SNR_DB/10);
  W = inv(H'*H + 1/sigm * eye(2))*H';
  mmse_rcv_beam = W * mimo_rcv;
  % signal detect
  mmse_rcv_beam = real(mmse_rcv_beam);
  % rcv_beam = rcv_beam / max(rcv_beam,[],'all');
  % mmse_rcv_beam(rcv_beam>0.3) = 1;
  % mmse_rcv_beam(rcv_beam<-0.3) = -1;
  % symbol error rate
  % mmse_symbol_err = 1-mean(sig==rcv_beam, 'all')
  mmse_snr = SNR(sig,mmse_rcv_beam);
end

function solution = add_awgn(X,SNR)
    N=length(X);
    signalPower = sum(X(1:end).^2)/N;
    linearSNR = 10^(SNR/10);
    [a,b]=size(X);
    noiseMat = randn(a,b).*sqrt(signalPower/linearSNR);
    solution = X + noiseMat;
end

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
