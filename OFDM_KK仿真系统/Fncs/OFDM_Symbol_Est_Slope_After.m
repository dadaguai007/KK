
% 接收信号
rxTrainSymbol_phase = data_kk(1:nTrainCarrier,:);
%参考信号
refTrainSymbol_phase = qam_signal_mat(1:nTrainCarrier,:);
% 计算角度
phi2=angle(rxTrainSymbol_phase./refTrainSymbol_phase);
% 对符号的相位偏差估计
S_est=zeros(size(qam_signal_mat,2),1);
for j=1:size(qam_signal_mat,2)
    K=zeros((nTrainCarrier),1);
    M=zeros((nTrainCarrier),1);
    for i=1:(nTrainCarrier)
        K(i)=i*phi2(i,j);
        M(i)=i.^2;
    end
    S_est(j)=sum(K)./sum(M);
end

figure;
plot(S_est)
xlabel('Index of OFDM symbol')