% 对每个符号下的所有子载波进行补偿

for k=1:size(qam_signal_mat,2)
    for m=1:size(qam_signal_mat,1)
        data_kk(m,k)=data_kk(m,k)*exp(-1j*m*S_l(k));
    end
end
