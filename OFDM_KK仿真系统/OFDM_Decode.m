%% kk
nTrainSym =100*k;% 频域均衡的参数
rxsig = ipd_btb(1:2*floor(length(ipd_btb)/2)).';

% 实施KK
c=0;
fup=fs*6;
[SSB_Sig_KK,ln_sig] = KK_MPC_Multi(rxsig+c,fs,fup);

%下采样
data_kk = downsample(SSB_Sig_KK,fup/fs);
% DATA=data_kk;
% % 对齐，滤波操作
% dataout1=NTF(DATA,fs,5e3,500e3,0);
% data_kk=dataout1;


% 解OFDM
data_kk_ofdm = reshape(data_kk,nn.fft_size+nn.nCP,[]);
data_kk_ofdm(1:nn.nCP,:) = [];
data_kk_ofdm = fft(data_kk_ofdm);
% get the modulated carriers
data_kk_ofdm = data_kk_ofdm(postiveCarrierIndex,:);
% channel estimation
rxTrainSymbol = data_kk_ofdm(:,1:nTrainSym);
qam_signal_mat=repmat(qam_signal,1,k);
refTrainSymbol = qam_signal_mat(:,1:nTrainSym);
Hf = mean(rxTrainSymbol./refTrainSymbol,2);

% channel equalization
data_kk = data_kk_ofdm.*repmat(1./Hf,1,nn.nPkts*k);

% 计算符号相位偏差
if strcmp(Sym_EST,'symbol_est')
    OFDM_Symbol_Est_Slope;
end

% 相位偏差消除
if strcmp(P_EST,'phase')
    phase_compensation;
end

if strcmp(f_EST,'fre_est')
    % OFDM_Phase_com;
    OFDM_Symbol_Est_Slope_After;
end

%保留信号矩阵
data_kk_mat=data_kk;
%归一化
data_kk=data_kk(:);
data_kk = data_kk./sqrt(mean(abs(data_kk(:)).^2));



%% 解码，计算误码
% 参考序列
ref_seq_1 =qamdemod(ref_seq,nn.M,'OutputType','bit','UnitAveragePower',1);
ref_seq_1=ref_seq_1(:);
% 接收序列
yyy = data_kk;

yyy_1 = qamdemod(yyy,nn.M,'OutputType','bit','UnitAveragePower',1);
yyy_1=yyy_1(:);
% 全部的序列进行解码
[ber1(index),num1(index),error_location] = CalcBER(yyy_1,ref_seq_1); %计算误码率
fprintf('Num of Errors = %d, BER = %1.7f\n',num1(index),ber1(index));

% 按符号解码，每个符号上的错码数
Calc_BER_mat;