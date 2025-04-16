% dither 对 多载波信号的影响
clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Project\Base_Code\Base\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')

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


    % IQ调制

    % 输入光信号
    Pi_dBm = 10;
    Pi = 10^(Pi_dBm/10)*1e-3; %W

    phi=0.87; % Vbias 偏移程度
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
    Cal_CSPR(sigTxo,Ai,Vpi,phi);

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
    [ReceivedSignal,~]=Receiver.Total_Preprocessed_signal(ipd_btb);
    % 归一化
    ReceivedSignal=pnorm(ReceivedSignal);
    Dc=mean(ReceivedSignal);
    % BER 计算
    [ber_total,num_total]=Receiver.Cal_BER(ReceivedSignal);


end

% KK恢复后频谱
[ IfdBm_Dith_dsb_receiver, fre ]=mon_ESA_flag(real(ReceivedSignal),fs,0);
% PD接收后电流
[ IfdBm_Dith_dsb_I, fre ]=mon_ESA_flag(ipd_btb,fs,0);

% 两者相减
IfdBm_Dith_dsb_II=(IfdBm_Dith_dsb_I)-(IfdBm_Dith_dsb_receiver);

% PD接收后电流
% [ IfdBm_Dith_dsb_II, fre ]=mon_ESA_flag(II,fs,0);


WB = OCG_WaitBar(k);
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
    [receive,Dc]=Receiver.Preprocessed_signal(ipd_btb);
    [signal_ofdm_martix,data_ofdm_martix,Hf,data_qam,qam_bit]=Receiver.Demodulation(receive);
    % BER 计算
    [ber,num]=Receiver.Cal_BER(receive);
    WB.updata(i);
end
WB.closeWaitBar();% 分段解码


% 选取性能较好段，进行重新调制
[ofdm_signal,~] = nn.ofdm(data_ofdm_martix);
% 补上直流
Re_Signal=Dc+ofdm_signal;
% 信号复制
signal_Re=repmat(Re_Signal,k,1);
% 残留噪声
Re=ReceivedSignal-signal_Re.';
% 更正解码参数
Receiver.Total_Preprocessed_signal(ipd_btb);
% 解码(接收信号减去残留噪声)
[ber_total1,num_total1]=Receiver.Cal_BER(ReceivedSignal-Re);

% 噪声功率谱

[ IfdBm, Fre ]=mon_ESA_flag(Re,fs,0);
figure;
plot(Fre,IfdBm,'b');
FontSize=14;
flag=struct();
flag.LegendON_OFF=0;
xticks([-32,-20,-10,0,10,20,32]);
Plotter('Residual noise spectrum','Frequency (GHz)','Magnitude (dB)',[-fs/2/1e9 fs/2/1e9],...
    [-150 50],'',flag,FontSize);

% 创建时间轴
[~,t_up]=freq_time_set(length(signal),fs);

figure;
plot(t_up,real(ReceivedSignal-Re),'b')
flag=struct();
flag.LegendON_OFF=0;
xticks([0,1e-5,2e-5,3e-5,3.5e-5]);
Plotter('Recovered signal','Time','Amplitude',[0 3.5e-5],[-1.5 1.5],...
    '',flag,FontSize);



% KK恢复后的信号谱
color=distinguishable_colors(20);
figure;hold on;
plot(fre,IfdBm_Dith_dsb_receiver,'Color',color(12,:));
plot(fre,IfdBm_Dith_dsb_I,'Color',color(17,:));

xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('DSB - After KK / After PD')
legend('KK','PD');
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 14);
set(gcf,'Position', [0, 0, 480, 400]);
set(gca, 'LineWidth', 1.25);
set(gca,'XLim',[-fs/2/1e9 fs/2/1e9],'YLim',[-150 30]);