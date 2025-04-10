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
for index= 1:length(flag_mon)
    % dither 信号 选择不同的dither方案
    type=char(flag_mon(index));
    % 方案选择
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

    Signal_Dither_Beat(index,:)=S+S1;
    % 频谱存储
    [ IfdBm, Fre ]=mon_ESA_flag(Signal_Dither_Beat(index,:),fs,0);
    PowerdBm(index,:)=IfdBm;
end



% 噪声功率谱

figure;hold on;
plot(Fre,PowerdBm(1,:),'b');
plot(Fre,PowerdBm(2,:),'r');
FontSize=14;
flag=struct();
flag.LegendON_OFF=1;
flag.Legendflage=1;
% 标题
LegendArray=["DSB","SSB"];
xticks([-32,-20,-10,0,10,20,32]);
Plotter('Signal-dither beat noise spectrum','Frequency (GHz)','Magnitude (dB)',[-fs/2/1e9 fs/2/1e9],...
    [-150 50],LegendArray,flag,FontSize);
