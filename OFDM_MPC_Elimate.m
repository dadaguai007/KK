% dither 对 多载波信号的影响
clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Codebase\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
% addpath('OFDM_KK仿真系统\Fncs\')
% addpath('OFDM_KK仿真系统\')
OFDM_TX;
% 生成信号
[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);
SpS=nn.Fs/nn.Rs;
% 采样率
fs=nn.Fs;

% 归一化
scale_factor = max(max(abs(real(signal))),max(abs(imag(signal))));
signal = signal./scale_factor;

%参考信号
ref_seq=reshape(qam_signal,1,[]);
ref_seq = repmat(ref_seq,1,100);
% 重复信号
k=20;
% qam信号矩阵
ref_seq_mat=repmat(qam_signal,1,k);

% 信号复制
signal=repmat(signal,k,1);
% dither 的频率处理
f1=40e3;
f2=60e3;
Fs_new=nn.Fs;
N=length(signal)/(Fs_new/f1);

%方案选择
type='dsb';


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


%IQ
%信号幅度
As=1;
% Dither 的幅度设置
alfa1=0.2;
alfa2=alfa1;
V1=alfa1*As;% I路dither幅度 双边带下两路的dither幅度
V2=alfa2*As;% Q路dither幅度
V=[V1,V2];

% vpp=num2str(alfa1);% 百分之Vpi


%IQ
N=length(label);

% 转置
signal=signal.';

%  全部解码
if 1
    Vdither=V;
    if strcmp(type,'dsb')
        VbI=Vdither(1)*Vdither1;
        VbQ=Vdither(1)*Vdither2;
    elseif strcmp(type,'ssb')
        VbI=Vdither(1)*Vdither1+Vdither(2)*Vdither_1;
        VbQ=Vdither(1)*Vdither2+Vdither(2)*Vdither_2;
    end

    % 载波幅度
    Ac=1;
    % 加载信号
    sigTxo=real(signal)+VbI+1j*imag(signal)+1j*VbQ;
    signal_power=signalpower(sigTxo);
    fprintf('optical signal power: %.2f dBm\n', 10 * log10(signal_power / 1e-3));

    % 传输信号
    sigTxo_DC=Ac+real(signal)+VbI+1j*imag(signal)+1j*VbQ;
    % dither
    Dither=VbI+1j*VbQ;
    % Cal CSPR
    Carry_power=Ac.^2;
    CSPR = Carry_power/signal_power;
    CSPR=10*log10(CSPR);
    CSPR=round(CSPR,2);
    text='CSPR=';
    leg_text=strcat(text,num2str(CSPR),'dB');
    fprintf("the CSPR is %1.2fdB\n",CSPR);

    %接收
    sigRxo=sigTxo_DC;

    %PD
    paramPD.B =nn.Fs;
    paramPD.R =1;
    paramPD.type = 'ideal';
    paramPD.Fs=nn.Fs;
    %pd
    ipd_btb = pd(sigRxo, paramPD);
    % dither 平方项
    ipd_pingfang = pd(Dither, paramPD);
    % dither的扩大
    S2=Dither*Ac+Ac*conj(Dither);
    % 拍频项 1
    S=As*signal.*VbI+conj(As*signal).*VbI;
    % 拍频项 2
    S1=real(As*signal.*conj(1j*VbQ)+conj(As*signal).*(1j*VbQ));
    % 去除干扰项
    %     ipd_btb=ipd_btb-S-S1-S2;
    % 发射机参数
    ofdmPHY=nn;
    %%---------------------------------------        解码       ---------------------------%%
    Receiver=OFDMreceiver( ...
        ofdmPHY, ...       %%% 发射机传输的参数
        ofdmPHY.Fs, ...    %   采样
        6*ofdmPHY.Fs, ...  % 上采样
        10*k, ...            % 信道训练长度
        1:1:ofdmPHY.nModCarriers, ...    %导频位置
        1, ...             % 选取第一段信号
        ref_seq, ...       % 参考序列
        ref_seq_mat, ...    % qam 矩阵
        'off', ...         % 是否采用CPE
        'off', ...         % 对所有载波进行相位补偿
        'KK');             % 接收方式

    % 信号预处理
    [ReceivedSignal,Dc]=Receiver.Total_Preprocessed_signal(ipd_btb);
    %     [signal_ofdm_martix,data_ofdm_martix,Hf,data_qam,qam_bit]=Receiver.Demodulation(ReceivedSignal);
    % BER 计算
    [ber_total,num_total]=Receiver.Cal_BER(ReceivedSignal);


end

mon_ESA(ReceivedSignal,fs);


WB = OCG_WaitBar(k);
% 发射机参数
ofdmPHY=nn;
for i=1:k

    %%---------------------------------------        解码       ---------------------------%%
    Receiver=OFDMreceiver( ...
        ofdmPHY, ...       %%% 发射机传输的参数
        ofdmPHY.Fs, ...    %   采样
        6*ofdmPHY.Fs, ...  % 上采样
        10, ...            % 信道训练长度
        1:1:ofdmPHY.nModCarriers, ...    %导频位置
        i, ...             % 选取第一段信号
        ref_seq, ...       % 参考序列
        qam_signal, ...    % qam 矩阵
        'off', ...         % 是否采用CPE
        'off', ...         % 对所有载波进行相位补偿
        'KK');             % 接收方式

    % 信号预处理
    [receive,Dc]=Receiver.Preprocessed_signal(ipd_btb);
    [signal_ofdm_martix,data_ofdm_martix,Hf,data_qam,qam_bit]=Receiver.Demodulation(receive);
    % BER 计算
    [ber,num]=Receiver.Cal_BER(receive);
    WB.updata(i);
end
WB.closeWaitBar();% 分段解码

% 选取性能较好段，进行重新调制
[ofdm_signal,~] = nn.ofdm(data_ofdm_martix);
Re_Signal=Ac+label;
% 信号复制
signal_Re=repmat(Re_Signal,k,1);
% 残留噪声
Re=ReceivedSignal-signal_Re.';
% 更正解码参数
Receiver.Total_Preprocessed_signal(ipd_btb);
% 解码
[ber_total1,num_total1]=Receiver.Cal_BER(ReceivedSignal-Re);

% 噪声功率谱
mon_ESA(Re,fs);