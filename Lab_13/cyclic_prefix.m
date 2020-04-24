% this function adds the last 36 bits of the symbole to the begining as cyclic prefix
function c_sig = cyclic_prefix(sig,cp_num)
  L = length(sig)
  c_sig = [sig(L-cp_num+1:end);sig]
end
