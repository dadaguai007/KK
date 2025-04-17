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

% %方案选择
% type='dsb';



%IQ
%信号幅度
As=1;
% Dither 的幅度设置
alfa1=0.2;
alfa2=alfa1;
V1=alfa1*As;% I路dither幅度 双边带下两路的dither幅度
V2=alfa2*As;% Q路dither幅度
V=[V1,V2];


% 转置
signal=signal.';

flag_mon={'dsb','ssb'};
WB = OCG_WaitBar(2);
%  全部解码
for index= 1:length(flag_mon)
    % dither 信号 选择不同的dither方案
    type=char(flag_mon(index));

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

    % 幅度
    Vdither=V;
    if strcmp(type,'dsb')
        VbI=Vdither(1)*Vdither1;
        VbQ=Vdither(1)*Vdither2;
    elseif strcmp(type,'ssb')
        VbI=Vdither(1)*Vdither1+Vdither(2)*Vdither_1;
        VbQ=Vdither(1)*Vdither2+Vdither(2)*Vdither_2;
    end

    % 载波幅度
    Ac=1.15;
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
%     text='CSPR=';
%     leg_text=strcat(text,num2str(CSPR),'dB');
    fprintf("the CSPR is %1.2fdB\n",CSPR);



    %无串扰的传输信号
    sigTxo_Pure=Ac+(signal);

    signal_power1=signalpower(signal);
    % Cal CSPR
    Carry_power1=Ac.^2;
    CSPR1 = Carry_power1/signal_power1;
    CSPR1=10*log10(CSPR1);
    CSPR1=round(CSPR1,2);
%     text='CSPR=';
%     leg_text=strcat(text,num2str(CSPR1),'dB');
    fprintf("the Pure CSPR is %1.2fdB\n",CSPR1);
   
    
    %接收
    sigRxo=sigTxo_DC;



    %PD
    paramPD.B =nn.Fs;
    paramPD.R =1;
    paramPD.type = 'ideal';
    paramPD.Fs=nn.Fs;
    %pd
    ipd_btb = pd(sigRxo, paramPD);

    % 无串扰接收
    ipd_btb_Pure = pd(sigTxo_Pure, paramPD);


    % dither 平方项
    ipd_pingfang = pd(Dither, paramPD);
    % dither的扩大
    S2=Dither*Ac+Ac*conj(Dither);
    % 拍频项 1
    S=As*signal.*VbI+conj(As*signal).*VbI;
    % 拍频项 2
    S1=real(As*signal.*conj(1j*VbQ)+conj(As*signal).*(1j*VbQ));

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
    ReceivedSignal=pnorm(ReceivedSignal);
    % BER 计算
    [ber_total,num_total]=Receiver.Cal_BER(ReceivedSignal);




    % 无串扰信号预处理
    [ReceivedSignal_Pure,Dc1]=Receiver.Total_Preprocessed_signal(ipd_btb_Pure);
    ReceivedSignal_Pure=pnorm(ReceivedSignal_Pure);
    % BER 计算
    [ber_total1,num_total1]=Receiver.Cal_BER(ReceivedSignal_Pure);

    if strcmp(type,'ssb')
        % 变量名

        [ IfdBm_Dith_ssb_receiver, fre ]=mon_ESA_flag(ReceivedSignal,fs,0);
    elseif strcmp(type,'dsb')

        [ IfdBm_Dith_dsb_receiver, fre ]=mon_ESA_flag(ReceivedSignal,fs,0);

    end
WB.updata(index);
end

WB.closeWaitBar();


[ IfdBm__receiver, fre ]=mon_ESA_flag(ReceivedSignal_Pure,fs,0);





% 创建时间轴
[~,t_up]=freq_time_set(length(signal),fs);
% 接收信号
figure;
plot(t_up,ipd_btb)
FontSize=14;
flag=struct();
flag.LegendON_OFF=0;
Plotter('Remove signal-dither and carrier-dither beat','Time','Amplitude',[0 3.5e-5],[-1 6],'',flag,FontSize);

% 恢复信号
figure;
plot(t_up,real(ReceivedSignal))
FontSize=14;
flag=struct();
flag.LegendON_OFF=0;
Plotter('After KK relation','Time','Amplitude',[0 3.5e-5],[-1 3],'',flag,FontSize);





% dither信号-无串扰信号频率对比
color=distinguishable_colors(20);
figure;hold on;
plot(fre,IfdBm__receiver,'Color',color(12,:));
plot(fre,IfdBm_Dith_dsb_receiver,'Color',color(2,:));
plot(fre,IfdBm_Dith_ssb_receiver,'Color',color(6,:));
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('After KK')
legend('Wo Dither','DSB','SSB');
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 14);
set(gcf,'Position', [0, 0, 480, 400]);
set(gca, 'LineWidth', 1.25);
set(gca,'XLim',[-fs/2/1e9 fs/2/1e9],'YLim',[-200 30]);