clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Project\Base_Code\Base\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
addpath("Plot\")
addpath('GUI\')

OFDM_TX;
% 生成信号
[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);
SpS=nn.Fs/nn.Rs;
% 采样率
fs=nn.Fs;

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
f1=400e3;
f2=600e3;
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
Vpi=paramIQ.Vpi;
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
    100, ...            % 信道训练长度
    1:1:ofdmPHY.nModCarriers, ...    %导频位置
    1, ...             % 选取第一段信号
    ref_seq, ...       % 参考序列
    ref_seq_mat, ...    % qam 矩阵
    'off', ...         % 是否采用CPE
    'off', ...         % 对所有载波进行相位补偿
    'KK', ...          % 接收方式
    'on');             % 是否全部接收

% 初始化设置
Eb_N0_dB=15:30;
ber_total=zeros(length(Eb_N0_dB),1);
num_total=zeros(length(Eb_N0_dB),1);
WB = OCG_WaitBar(length(Eb_N0_dB));
for index=1:length(Eb_N0_dB)
    fprintf('噪声为 %d dB\n',Eb_N0_dB(index))
    % 噪声
    noise=EbN0_dB(sigRxo,Eb_N0_dB(index));

    % 加入噪声
    pd_receiver = real(pnorm(ipd_btb)+noise);

    % 对信号进行切分，并提出全部信号
    [DataGroup,totalPortion]=Receiver.Synchronization(pd_receiver);

    % 信号预处理
    [ReceivedSignal,dc]=Receiver.Preprocessed_signal(totalPortion);
    Dc=mean(ReceivedSignal);
    % BER 计算
    [ber_total(index),num_total(index)]=Receiver.Cal_BER(ReceivedSignal);

    % 分组kk
    re_signal=[];
    for Idx=1:k
        % 序列号
        Receiver.Nr.squ_num=Idx;
        % k设置为1
        Receiver.Nr.k=1;
        selectedPortion=Receiver.selectSignal(k,DataGroup);
        % KK
        [selectSignal,dc]=Receiver.Preprocessed_signal(selectedPortion);
        % BER 计算
        [ber(Idx),num(Idx),l(Idx)]=Receiver.Cal_BER(selectSignal);
        % 存储kk
        re_signal=[re_signal,selectSignal];
    end

    fprintf('分组解码的BER = %1.7f\n',sum(num)/sum(l));
    ber_group_total1(index)=sum(num)/sum(l);



    if 1
        disp('进入循环');

        % 更新PD输入
        pd_input=pd_receiver;
        for jj=1:20
            % 算法迭代
            alpha=0.02;
            recoverI=iteraElimate(re_signal,pd_input,fs,alpha,real(mean(re_signal)),Vdither(1));
            % 对信号进行切分，并提出全部信号
            [DataGroup1,totalPortion1]=Receiver.Synchronization(recoverI);
            % 分组kk
            re_signal1=[];
            for Idx=1:k
                % 序列号
                Receiver.Nr.squ_num=Idx;
                % k设置为1
%                 Receiver.Nr.k=1;
                selectedPortion1=Receiver.selectSignal(k,DataGroup1);
                % KK
                [selectSignal1,~]=Receiver.Preprocessed_signal(selectedPortion1);
                % 存储kk
                re_signal1=[re_signal1,selectSignal1];
            end
            % 赋值
            re_signal=re_signal1;
            % BER 计算
            [ber_total1(jj),num_total1]=Receiver.Cal_BER(re_signal);
            % 更新PD输入
            pd_input=recoverI;
            martixI(jj,:)=recoverI;
        end
        [ber_total_iter(index),inn]=min(ber_total1);

        % 进行拍频的迭代


        % 更新PD输入
        pd_input=martixI(inn,:);
        for jj=1:20
            % 算法迭代
            alpha=0.02;
            recoverI=iteraElimate_beat(re_signal,pd_input,fs,alpha,real(mean(re_signal)),Vdither(1));
            % 对信号进行切分，并提出全部信号
            [DataGroup1,totalPortion1]=Receiver.Synchronization(recoverI);
            % 分组kk
            re_signal1=[];
            for Idx=1:k
                % 序列号
                Receiver.Nr.squ_num=Idx;
                % k设置为1
%                 Receiver.Nr.k=1;
                selectedPortion1=Receiver.selectSignal(k,DataGroup1);
                % KK
                [selectSignal1,~]=Receiver.Preprocessed_signal(selectedPortion1);
                % 存储kk
                re_signal1=[re_signal1,selectSignal1];
            end
            % 赋值
            re_signal=re_signal1;
            % BER 计算
            [ber_total2(jj),num_total1]=Receiver.Cal_BER(re_signal);
            % 更新PD输入
            pd_input=recoverI;
           
        end
        ber_total_iter2(index)=min(ber_total2);

    end
    WB.updata(index);

end
WB.closeWaitBar();



berplot = BERPlot_David();
% 间隔
berplot.interval=3;
% 字号
berplot.Config.FontSize = 14;
berplot.flagThreshold=1;
berplot.flagRedraw=0;
berplot.flagAddLegend=1;
% BER=[ber_total.';ber_total_iter;ber_total_iter_Group];
BER=[ber_total.';ber_group_total1;ber_total_iter2];
LengendArrary=["80km w/o ","80km w/o Group",'80km w SIC'];
berplot.multiplot(Eb_N0_dB,BER,LengendArrary);
