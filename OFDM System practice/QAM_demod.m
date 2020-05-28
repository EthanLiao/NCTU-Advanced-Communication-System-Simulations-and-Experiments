function demod_sig = QAM_demod(sig)
  x = (0:15)
  q_x = QAM(x)

  for i=1:length(sig)
    [minVal,minArg] =  min(abs(sig(i)-q_x))
    demod_sig(i) = x(minArg)
end
