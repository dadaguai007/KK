clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Codebase\')
addpath('OFDM_KK仿真系统\Fncs\')
addpath('OFDM_KK仿真系统\')
OFDM_TX;
% 生成信号
[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);
SpS=nn.Fs/nn.Rs;
% 采样率
fs=nn.Fs;


scale_factor = max(max(abs(real(signal))),max(abs(imag(signal))));
signal = signal./scale_factor;

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
type='dsb';
pre='_ssb';
flagtag=0;% 绘BER图标致
savetag=0;% 存储标志
savefigtag=0;% 存图标志
showtag='off';% 图片是否显示


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
%信号幅度
As=1;
% Dither 的幅度设置
alfa1=0.1;
alfa2=alfa1;
V1=alfa1*As;% I路dither幅度 双边带下两路的dither幅度
V2=alfa2*As;% Q路dither幅度
V=[V1,V2];

vpp=num2str(alfa1);% 百分之Vpi

if strcmp(Amp_type,'sweep_Amp')
    afla=0.01:0.01:0.1; % 比例
    for hindex=1:length(afla)
        for hindex1=1:length(afla)
            V2(hindex1)=afla(hindex1)*9;
            v1(hindex,hindex1)=afla(hindex)*9;
        end
    end
    % 转置为列
    v1=v1.';
    V1=v1(:);
    V2=V2.';
    V2=repmat(V2,10,1);
    V=[V1,V2];
end


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

signal=signal.';

snr=25;
WB = OCG_WaitBar(length(snr));
for Index=1:size(V,1)
    Vdither=V(Index,:);

    for index=1:length(snr)

        if strcmp(type,'dsb')
            VbI=Vdither(1)*Vdither1;
            VbQ=Vdither(1)*Vdither2;
        elseif strcmp(type,'ssb')
            VbI=Vdither(1)*Vdither1+Vdither(2)*Vdither_1;
            VbQ=Vdither(1)*Vdither2+Vdither(2)*Vdither_2;
        end



        % 载波幅度
        Ac=0.9;
        % 加载信号
        sigTxo=real(signal)+VbI+1j*imag(signal)+1j*VbQ;
        signal_power=signalpower(sigTxo);
        fprintf('optical signal power: %.2f dBm\n', 10 * log10(signal_power / 1e-3));

        sigTxo_DC=Ac+real(signal)+VbI+1j*imag(signal)+1j*VbQ;

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

        OFDM_Decode;
        OFDM_EVM;
        if strcmp(Amp_type,'sweep_Amp')
            BER(Index,index)=ber1(index);
            CSPR_mat(Index,index)=CSPR;
        else
            BER(index)=ber1(index);
            CSPR_mat(index)=CSPR;
        end
        WB.updata(index);
    end
end
WB.closeWaitBar();





if savetag
    Title='CSPR_out';
    datapath = strcat('Output\',Title,pre,vpp);
    if ~exist(datapath,'dir')
        mkdir(datapath);
    end
    save(sprintf('%s\\CSPR_mat.mat',datapath),'CSPR_mat');


    Title='Averagepower_out';
    datapath = strcat('Output\',Title,pre,vpp);
    if ~exist(datapath,'dir')
        mkdir(datapath);
    end
    save(sprintf('%s\\Average_mat.mat',datapath),'averagepower_mat');


    Title='BER_out';
    datapath = strcat('Output\',Title,pre,vpp);
    if ~exist(datapath,'dir')
        mkdir(datapath);
    end
    save(sprintf('%s\\BER_mat.mat',datapath),'BER');
end

if flagtag
    berplot = BERPlot_David();

    berplot.interval=2;
    berplot.flagThreshold=1;
    berplot.flagRedraw=1;
    berplot.flagAddLegend=0;
    berplot.orderInterpolation=2;
    berplot.plot(snr,BER);
end