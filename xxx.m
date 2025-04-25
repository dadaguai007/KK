% 测试dither  与 信号的 损伤比

clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Project\Base_Code\Base\')
% addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
addpath('Plot\')
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
As=9;
% Dither 的幅度设置
alfa1=0.02;
alfa2=alfa1;
V1=alfa1*As;% I路dither幅度 双边带下两路的dither幅度
V2=alfa2*As;% Q路dither幅度
V=[V1,V2];


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


    % IQ调制

    % 输入光信号
    Pi_dBm = 10;
    Pi = 10^(Pi_dBm/10)*1e-3; %W

    phi=0.86; % Vbias 偏移程度
    %     phi=0.83; % 对应0.05
    paramIQ.Vpi=9;
    if strcmp(type,'dsb')
        paramIQ.VbI=-phi*paramIQ.Vpi+VbI;
        paramIQ.VbQ=-phi*paramIQ.Vpi+VbQ;
    elseif strcmp(type,'ssb')
        paramIQ.VbI=-phi*paramIQ.Vpi+VbI;
        paramIQ.VbQ=-phi*paramIQ.Vpi+VbQ;
    end
    paramIQ.Vphi=paramIQ.Vpi/2;
    Amp=1.7; % 信号放大
    Ai  = sqrt(Pi);% 输入光功率
    m=Modulation_index(Amp*signal.',paramIQ.Vpi,'ofdm');
    fprintf('Modulation index=%3.3f\n',m);
    % 调制
    sigTxo = iqm(Ai, Amp*signal, paramIQ);


    signal_power=signalpower(sigTxo);
    fprintf('optical signal power: %.2f dBm\n', 10 * log10(signal_power / 1e-3));

    % Cal CSPR
    Vpi=9;
    Eout1=Cal_CSPR(sigTxo,Ai,Vpi,phi);


    V=VbI+1j*VbQ;
    I=signalpower(V);

    P=Amp*signal.';
    Q=signalpower(P);
    pental = I/Q;
    pental_dB=10*log10(pental);

    fprintf('dither to signal pental power: %.2f dB\n', pental_dB);
    % 光纤传输
    % fiber param
    param=struct();
    param.Ltotal = 80; %km
    param.Lspan =10;
    param.hz= 0.5;
    param.alpha=0.2;
    param.D = 16;
    param.gamma = 1.3;
    param.Fc = 193.1e12;
    param.NF = 4.5;
    param.amp='ideal';
    param.Fs=fs;

    % OFDM_tran
    sigRxo=ssfm(sigTxo,param);
    power2=signalpower(sigRxo);
    fprintf('after ssfm signal power: %.2f dBm\n', 10 * log10(power2 / 1e-3));


    %PD
    paramPD.B =nn.Fs;
    paramPD.R =1;
    paramPD.type = 'ideal';
    paramPD.Fs=nn.Fs;
    %pd
    ipd_btb = pd(sigRxo, paramPD);
    ipd_btb=pnorm(ipd_btb);

    Eb_N0_dB=15:30;
    for index=1:length(Eb_N0_dB)
        % 噪声
        noise=EbN0_dB(sigRxo,Eb_N0_dB(index));

        % 加入噪声
        pd_receiver = real(ipd_btb+noise);

        re_signal=[];
        % 发射机参数
        ofdmPHY=nn;
        for i=1:k

            %%---------------------------------------        解码       ---------------------------%%
            Receiver=OFDMreceiver( ...
                ofdmPHY, ...       %%% 发射机传输的参数
                ofdmPHY.Fs, ...    %   采样
                6*ofdmPHY.Fs, ...  % 上采样
                ofdmPHY.nPkts, ...            % 信道训练长度
                1:1:ofdmPHY.nModCarriers, ...    %导频位置
                i, ...             % 选取第一段信号
                ref_seq, ...       % 参考序列
                qam_signal, ...    % qam 矩阵
                'off', ...         % 是否采用CPE
                'off', ...         % 对所有载波进行相位补偿
                'KK');             % 接收方式

            % 信号预处理
            [receive,~]=Receiver.Preprocessed_signal(pd_receiver);
            % BER 计算
            [ber,num(i),l(i)]=Receiver.Cal_BER(receive);
            re_signal=[re_signal,receive];

        end

        fprintf('分组解码的BER = %1.7f\n',sum(num)/sum(l));
        ber_group(index)=sum(num)/sum(l);
        % 用分组出来的信号代替
        re=re_signal;
        Dc=mean(re_signal);
    end

end
load('Output\BER_Algriom_Dither_decode_one_by_one_0_2_5_8_10.mat')
load('Output\BER_Algriom_Dither_decode_Group_0_2_5_8_10.mat')
load('Output\BER_Dither_0_2_5_8_10.mat')
berplot = BERPlot_David();
% 间隔
berplot.interval=3;
% 字号
berplot.Config.FontSize = 14;
berplot.flagThreshold=1;
berplot.flagRedraw=0;
berplot.flagAddLegend=1;
% BER=[ber_total.';ber_total_iter;ber_total_iter_Group];
BER=[ber_total(1,:);...
    ber_group;ber_total_iter_Group(2,:);];
LengendArrary=["80km w/o 0%Vpi",...
    "80km w/o Group 2%Vpi","80km w SIC Group 2%Vpi",];
berplot.multiplot(Eb_N0_dB,BER,LengendArrary);