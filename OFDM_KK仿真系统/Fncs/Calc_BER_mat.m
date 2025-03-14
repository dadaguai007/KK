% OFDM 矩阵解码
ref_seq_mat1 =qamdemod(ref_seq_mat,nn.M,'OutputType','bit','UnitAveragePower',1);

% 接收序列

yyy_mat = qamdemod(data_kk_mat,nn.M,'OutputType','bit','UnitAveragePower',1);
% 按照OFDM符号解码,不能用作误码率计算
for idx=1:size(data_kk_mat,2)
    ref_symbol=ref_seq_mat1(:,idx);
    yyy_symbol=yyy_mat(:,idx);
    [ber,num,error_location] = CalcBER(yyy_symbol,ref_symbol);
    number_num(idx)=num;

end
