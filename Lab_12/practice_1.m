clear all;clf;close all;


N = 8;
signal = [0 1 2 3 4 5 6 7];
mod_signal = Npsk(signal,N);
scatterplot(mod_signal);

function mod_sig = Npsk(sig,N)
  for i = 1:length(sig)
      mod_sig(i) = cos(sig(i)*2*pi/N)+j*sin(sig(i)*2*pi/N);
  end
end
