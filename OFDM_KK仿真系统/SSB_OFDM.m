clear;close all;clc;
addpath('Fncs\')
addpath('D:\001-处理中\相干算法\optical_communication')
OFDM_TX;
% 生成信号
[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);
SpS=nn.Fs/nn.Rs;

fs=nn.Fs;
scale_factor = max(max(abs(real(signal))),max(abs(imag(signal))));
signal = signal./scale_factor;

% 记录最原始的数据长度
N_index=length(signal);

%参考信号
ref_seq=reshape(qam_signal,1,[]);
ref_seq = repmat(ref_seq,1,100);

% 重复信号
k=40;


% qam信号矩阵
ref_seq_mat=repmat(qam_signal,1,k);


signal=repmat(signal,k,1);
% dither 的频率处理
f1=40e3;
f2=60e3;
Fs_new=nn.Fs;
N=length(signal)/(Fs_new/f1);

%方案选择
type='ssb'; % 一个周期为5e4
Train_type='btb';% 传输方案选择

if strcmp(type,'dsb')
    Vdither1 = Creat_dither(Fs_new,f1,N);
    Vdither2 = Creat_dither(Fs_new,f2,N*(f2/f1));
elseif strcmp(type,'ssb')
    % f1 and f2 的cos信号
    Vdither1 = Creat_dither(Fs_new,f1,N);
    Vdither_1=Creat_dither(Fs_new,f2,N*(f2/f1));
    % f1 and f2 的sin信号
    Vdither_2 = Creat_dither1(Fs_new,f2,N*(f2/f1));
    Vdither2=Creat_dither1(Fs_new,f1,N);
end

Amp_type='none'; % 选择单边带的dither强度扫描

%IQ
%RF and dither
Amp=1.7;


% Dither 的幅度设置
alfa1=5e-2;
alfa2=alfa1;
V1=alfa1*9;% I路dither幅度 双边带下两路的dither幅度
V2=alfa2*9;% Q路dither幅度
V=[V1,V2];

vpp=num2str(alfa1);% 百分之Vpi


%IQ
N=length(label);
% 相位均衡参数
W=nn.nModCarriers;
pilotIndex=1:1:W;
P_EST='none';
% 相位估计
nTrainCarrier=nn.nModCarriers;% 估计符号的相位
Sym_EST='none';%symbol_est
f_EST='none';% 检查相位补偿后的符号相位分布
if 0
    CR =2.3; % Clipping Ratio
    signal_orgin=signal;
    signal = clipping(signal.',CR);
    figure;
    subplot(2,1,1)
    plot(real(signal_orgin))
    subplot(2,1,2)
    plot(real(signal))
else
    signal=signal.';
end

phi=0.87; % Vbias 偏移程度
% snr=5:2:25;% 系统的SNR
snr=25;
% WB = OCG_WaitBar(length(snr));
for Index=1:size(V,1)
    Vdither=V(Index,:);

    for index=1:length(snr)

        Pi_dBm = 10;
        Pi = 10^(Pi_dBm/10)*1e-3; %W

        paramIQ.Vpi=9;
        if strcmp(type,'dsb')
            paramIQ.VbI=-phi*paramIQ.Vpi+Vdither(1)*Vdither1;
            paramIQ.VbQ=-phi*paramIQ.Vpi+Vdither(1)*Vdither2;
        elseif strcmp(type,'ssb')
            paramIQ.VbI=-phi*paramIQ.Vpi+Vdither(1)*Vdither1+Vdither(2)*Vdither_1;
            paramIQ.VbQ=-phi*paramIQ.Vpi+Vdither(1)*Vdither2+Vdither(2)*Vdither_2;
        end
        paramIQ.Vphi=paramIQ.Vpi/2;

        Ai  = sqrt(Pi);
        m=Modulation_index(Amp*signal.',paramIQ.Vpi,'ofdm');
        fprintf('Modulation index=%3.3f\n',m);
        sigTxo = iqm(Ai, Amp*signal, paramIQ);

        %noise
        %         sigTxo=awgn(sigTxo,snr(index),'measured');

        %power
        power=signalpower(sigTxo);
        fprintf('optical signal power: %.2f dBm\n', 10 * log10(power / 1e-3));

        % fiber
        param=struct();
        param.Ltotal = 80; %km
        param.Lspan =80;
        param.hz= 0.1;
        param.alpha=0.2;
        param.D = 16;
        param.gamma = 1.3;
        param.Fc = 193.1e12;
        param.NF = 4.5;
        param.amp='ideal';
        param.Fs=nn.Fs;

        if strcmp(Train_type,'ssfm')
            OFDM_Train;
        else
            sigRxo=sigTxo;
        end
        %PD
        paramPD.B =nn.Fs;
        paramPD.R =1;
        paramPD.type = 'ideal';
        paramPD.Fs=nn.Fs;
        %pd
        ipd_btb = pd(sigRxo, paramPD);
        % 转换为电信号后加入趋势去除步骤
        % 平滑滤波的方式淘汰
        %         trand=smooth(ipd_btb);
        %         ipd_btb=ipd_btb-trand.';
        % 时域窗口平滑
%         size_window=10;
%         sig_rwin = rolling_window_central(ipd_btb, size_window, 'true');
%         ipd_btb=ipd_btb-mean(sig_rwin,2).';
%         dc=mean(ipd_btb);

        % PD
        %         ipd_btb=NTF(ipd_btb,fs,5e3,500e3,0);
        % 高通
        %         ipd_btb=real(HPF1(ipd_btb,fs,1e6,0));
        %         ipd_btb=ipd_btb+dc;


        % 按顺序解码
        OFDM_Decode_squence;

        % 全部解码
        %         OFDM_Decode;
        %         OFDM_EVM;
        CSPR_Mea;
        %         if strcmp(Amp_type,'sweep_Amp')
        %             BER(Index,index)=ber1(index);
        %             CSPR_mat(Index,index)=CSPR;
        %             averagepower_mat(Index,index)=U;
        %         else
        %             BER(index)=ber1(index);
        %             CSPR_mat(index)=CSPR;
        %             averagepower_mat(index)=U;
        %         end
        %         Plot_show;
        %         WB.updata(index);
    end
end
% WB.closeWaitBar();

% 频谱绘图
% mon_ESA_index(sigTxo,fs);
% mon_ESA_index(DATA,fs);
% mon_ESA_index(dataout1,fs);

scatplot(real(signal_scatter(:,19)),imag(signal_scatter(:,19)));

figure;
plot(symbol_EVM(:,19))
% figure;
% plot(subcarrier_EVM(:,20))
figure;
plot(subcarrier_index_symbol_EVM(:,19))

figure;
subplot(2,1,1)
plot(real(signal_scatter(1:2e5,19)),'.')

subplot(2,1,2)
plot(real(signal_squ(1:5e4,19)))