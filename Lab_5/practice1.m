clf;clear all;
% 4PAM Mapping
M_PAM = 8
x = (1:8)'
% y = pammod(x,M_PAM,pi/2)
y = PAMmod(x)
scatterplot(y)
text(real(y)+0.1,imag(y)-0.1,dec2bin(x))

% 16 QAM Mapping
M_QAM = 16  % modulation order
x = (0:15)' % generate a modulation symbol
y = QAMmod(x)
scatterplot(y)
text(real(y)+0.1,imag(y)+0.1,dec2bin(x))
title('16 QAM gray mapping')


% N = 3.0;
% x=linspace(-N, N);
% y=x;
% [X,Y]=meshgrid(x,y);
% z=(1000/sqrt(2*pi).*exp(-(X.^2/2)-(Y.^2/2)));
% surf(X,Y,z);
% shading interp
% axis tight

function c_sig = PAMmod(sig)
  for i=1:length(sig)
    if sig(i) == 1
      c_sig(i) = -7
      continue
    elseif sig(i) == 2
      c_sig(i) = -5
      continue
    elseif sig(i) == 3
      c_sig(i) = -3
      continue
    elseif sig(i) == 4
      c_sig(i) = -1
      continue
    elseif sig(i) == 5
      c_sig(i) = 1
      continue
    elseif sig(i) == 6
      c_sig(i) = 3
      continue
    elseif sig(i) == 7
      c_sig(i) = 5
      continue
    elseif sig(i) == 8
      c_sig(i) = 7
      continue
    end
  end
end

function c_sig = QAMmod(sig)
  sig = complex(sig)
  for i=1:length(sig)
    if sig(i)== complex(0)     % 0
      c_sig(i) = -3+3*j
      continue
    elseif sig(i)==complex(1) % 1
      c_sig(i) = -3+j
      continue
    elseif sig(i)==complex(2) % 2
      c_sig(i) = -3-3*j
      continue
    elseif sig(i)==complex(3) % 3
      c_sig(i) = -3-j
      continue
    elseif sig(i)==complex(4) % 4
      c_sig(i) = -1+3*j
      continue
    elseif sig(i)==complex(5) % 5
      c_sig(i) = -1+j
      continue
    elseif sig(i)==complex(6) % 6
      c_sig(i) = -1-3*j
      continue
    elseif sig(i)==complex(7) % 7
      c_sig(i) = -1-j
      continue
    elseif sig(i)==complex(8) % 8
      c_sig(i) = 3+3*j
      continue
    elseif sig(i)==complex(10) % 10
      c_sig(i) = 3-3*j
      continue
    elseif sig(i)==complex(9) % 10
      c_sig(i) = 3+j
      continue
    elseif sig(i)==complex(11) % 11
      c_sig(i) = 3-j
      continue
    elseif sig(i)==complex(12) % 12
      c_sig(i) = 1+3*j
      continue
    elseif sig(i)==complex(13) % 13
      c_sig(i) = 1+j
      continue
    elseif sig(i)==complex(14) % 14
      c_sig(i) = 1-3*j
      continue
    elseif sig(i)==complex(15) % 15
      c_sig(i) = 1-j
      continue
    end
  end
end
