clf;clear all;close all;
N = 1024;
sig = randi([0 1],2,N);
sig(sig==0) = -1;

h = randn(2,2);
H = h*(1/sqrt(2)+1/sqrt(2)*1j);

rcv_sig = H * sig;

SNR_DB = -20;
rcv_sig_1 = add_awgn_noise(rcv_sig(1,:), SNR_DB);
rcv_sig_2 = add_awgn_noise(rcv_sig(2,:), SNR_DB);
rcv_sig = [rcv_sig_1 ; rcv_sig_2];

% zero forcing
zf_rcv_sig = real(inv(H) * rcv_sig);
% zf_rcv_sig(zf_rcv_sig>0) = 1;
% zf_rcv_sig(zf_rcv_sig<0) = -1;
% zf_SER = 1-mean(zf_rcv_sig==sig, 'all')
zf_snr = SNR(sig, zf_rcv_sig)

sigm = 10^(SNR_DB/10);
% mmse
mmse_rcv_sig = real(inv(H'*H+1/sigm*eye(2))*H'*rcv_sig);
% mmse_rcv_sig(mmse_rcv_sig>0) = 1;
% mmse_rcv_sig(mmse_rcv_sig<0) = -1;
% mmse_SER = 1-mean(mmse_rcv_sig==sig, 'all')
mmse_snr = SNR(sig, mmse_rcv_sig)

subplot(2,2,1);stem(sig(1,:));title("mmse trans branch 1");
subplot(2,2,2);stem(sig(2,:));title("mmse trans branch 2");
subplot(2,2,3);stem(mmse_rcv_sig(2,:));title("mmse recieve branch 1");
subplot(2,2,4);stem(mmse_rcv_sig(2,:));title("mmse recieve branch 2");

figure();
subplot(2,2,1);stem(sig(1,:));title("zf trans branch 1");
subplot(2,2,2);stem(sig(2,:));title("zf trans branch 2");
subplot(2,2,3);stem(zf_rcv_sig(2,:));title("zf recieve branch 1");
subplot(2,2,4);stem(zf_rcv_sig(2,:));title("zf recieve branch 2");

function noise_sig = add_awgn_noise(sig, SNR_DB)
  L = length(sig);
  SYME = sum(abs(sig).^2)/L;
  SNR_n = 10^(SNR_DB/10);
  NYME = SYME / SNR_n;
  if isreal(sig)
    n = sqrt(NYME) .* randn(1,L);
  else
    n = sqrt(NYME/2) * (1+j) .* randn(1,L);
  end
  noise_sig = n+sig;
end

function snr = SNR(sig, rcv_sig)
  snr = 10*log10(mean(abs(sig).^2)/mean(abs(rcv_sig-sig).^2));
end
