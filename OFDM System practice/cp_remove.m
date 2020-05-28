function decp_data = cp_remove(sig,cp_sig_len,cp_num)
  recieve_mat = reshape(sig,cp_sig_len,ceil(length(sig)/cp_sig_len))
  decp_data = recieve_mat(cp_num+1:cp_sig_len)
end
