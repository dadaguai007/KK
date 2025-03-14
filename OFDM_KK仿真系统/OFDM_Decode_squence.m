% OFDM_Decode_squence
rxsig = ipd_btb(1:2*floor(length(ipd_btb)/2)).';
c=0;
fup=fs*4;
[SSB_Sig_KK,ln_sig] = KK_MPC_Multi(rxsig+c,fs,fup);


%下采样
data_kk = downsample(SSB_Sig_KK,fup/fs);
DATA=data_kk;
% DATA=sigRxo;
% 对齐，滤波操作
% dataout1=NTF(DATA,fs,5e3,500e3,0);
% DATA_filter=dataout1;

% dc_re=mean(DATA);
% dataout1=HPF1(DATA,fs,1e6,0);
% DATA_filter=dataout1+dc_re;
for squ_num=1:k
    % 训练序列
    nTrainSym=100;
    HK=1;
    % 选取每一个阶段的信号
    DATA_squ=DATA(N_index*(squ_num-1)+1:N_index*squ_num);

    %     % 去趋势
    %     trand=smooth(DATA_squ);
    %     DATA_squ=DATA_squ-trand.';

    % 滑动窗口函数，对数据进行排列
%             size_window=10;
%             sig_rwin = rolling_window(DATA_squ, size_window, 'true');
%             DATA_squ=DATA_squ.'-mean(sig_rwin,2);
    % 解OFDM
    data_kk_ofdm = reshape(DATA_squ,nn.fft_size+nn.nCP,[]);
    data_kk_ofdm(1:nn.nCP,:) = [];
    data_kk_ofdm = fft(data_kk_ofdm);
    % get the modulated carriers
    data_kk_ofdm = data_kk_ofdm(postiveCarrierIndex,:);
    % channel estimation
    rxTrainSymbol = data_kk_ofdm(:,1:nTrainSym);
    qam_signal_mat=repmat(qam_signal,1,HK);
    refTrainSymbol = qam_signal_mat(:,1:nTrainSym);
    Hf = mean(rxTrainSymbol./refTrainSymbol,2);

    % channel equalization
    data_kk = data_kk_ofdm.*repmat(1./Hf,1,nn.nPkts*HK);

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

    % 解码，计算误码
    % 参考序列
    ref_seq_1 =qamdemod(ref_seq,nn.M,'OutputType','bit','UnitAveragePower',1);
    ref_seq_1=ref_seq_1(:);
    % 接收序列
    yyy = data_kk;

    yyy_1 = qamdemod(yyy,nn.M,'OutputType','bit','UnitAveragePower',1);
    yyy_1=yyy_1(:);
    % 全部的序列进行解码
    [ber1(squ_num),num1(squ_num),error_location] = CalcBER(yyy_1,ref_seq_1); %计算误码率
    fprintf('Num of Errors = %d, BER = %1.7f\n',num1(squ_num),ber1(squ_num));
    % 按符号解码，每个符号上的错码数
    Calc_BER_mat;

    OFDM_EVM;
    symbol_EVM(:,squ_num)=rmsEVM_symbol.';
    subcarrier_EVM(:,squ_num)=rmsEVM_subcarrier;

    % 第L个符号的EVM
    LL=200;
    OFDM_symbol_EVM;
    subcarrier_index_symbol_EVM(:,squ_num)=rmsEVM_subcarrier_index_symbol;

    signal_scatter(:,squ_num)=data_kk;
    signal_squ(:,squ_num)=DATA_squ;

end
