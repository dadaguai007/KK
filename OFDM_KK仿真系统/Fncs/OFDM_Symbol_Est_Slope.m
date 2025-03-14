
% 接收信号
rxTrainSymbol_phase = data_kk(1:nTrainCarrier,:);
%参考信号
refTrainSymbol_phase = qam_signal_mat(1:nTrainCarrier,:);
% 计算角度
phi1=angle(rxTrainSymbol_phase./refTrainSymbol_phase);
% 对符号的相位偏差估计
S_est=zeros(size(qam_signal_mat,2),1);
for j=1:size(qam_signal_mat,2)
    K=zeros((nTrainCarrier),size(qam_signal_mat,2));
    M=zeros((nTrainCarrier),size(qam_signal_mat,2));
    for i=1:nTrainCarrier
        K(i,j)=i*phi1(i,j);
        M(i,j)=i.^2;
    end
    S_l(j)=sum(K(:,j))./sum(M(:,j));
end
% for v=1:x
%     for vv=1:size(qam_signal_mat,2)
%     K1(vv)=vv*S_l(vv);
%     K2(vv)=vv.^2;
%     end
%     Sn(v)=sum(K1)./sum(K2);
% end

    
% figure;
% plot(S_l)
% xlabel('Index of OFDM symbol')