% 找到最佳的alpha，通过光纤传输系统

clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Project\Base_Code\Base\')
% addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
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
Dc=mean(pnorm(sigTxo));
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
% 归一化
ipd_btb=pnorm(ipd_btb);



% 无dither信号
paramIQ.Vpi=9;
if strcmp(type,'dsb')
    paramIQ.VbI=-phi*paramIQ.Vpi;
    paramIQ.VbQ=-phi*paramIQ.Vpi;
elseif strcmp(type,'ssb')
    paramIQ.VbI=-phi*paramIQ.Vpi;
    paramIQ.VbQ=-phi*paramIQ.Vpi;
end
paramIQ.Vphi=paramIQ.Vpi/2;
Amp=1.7; % 信号放大
Ai  = sqrt(Pi);% 输入光功率
% 调制
sigTxo1 = iqm(Ai, Amp*signal, paramIQ);



% OFDM_tran
sigRxo1=ssfm(sigTxo1,param);
power3=signalpower(sigRxo1);
fprintf('after ssfm signal power: %.2f dBm\n', 10 * log10(power3 / 1e-3));







S=pd((sigRxo1),paramPD);
S=pnorm(S);

% 光电流差值
Error=ipd_btb-S;

figure;hold on;
plot(S)
plot(ipd_btb)
plot(Error)


% 误差绘图
figure
plot(Error)


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


pd_receiver = real(ipd_btb);

% 对信号进行切分，并提出全部信号
[DataGroup,totalPortion]=Receiver.Synchronization(pd_receiver);

% 信号预处理
[ReceivedSignal,dc]=Receiver.Preprocessed_signal(totalPortion);
Receiver.Button.Display='on';
% BER 计算
[ber_total,num_total]=Receiver.Cal_BER(ReceivedSignal);

% 分组kk
re_signal=[];
Receiver.Button.Display='off';
for Idx=1:k
    % 序列号
    Receiver.Nr.squ_num=Idx;
    % k设置为1
    Receiver.Nr.k=1;
    Receiver.Nr.nTrainSym=50;
    selectedPortion=Receiver.selectSignal(k,DataGroup);
    % KK
    [selectSignal,dc]=Receiver.Preprocessed_signal(selectedPortion);
    % 解码运算
    [signal_ofdm_martix,data_ofdm_martix,Hf,data_qam,qam_bit]=Receiver.Demodulation(selectedPortion);
    % 转换为矩阵形式
    data_mat=reshape(data_qam,nn.nModCarriers,[]);
    % BER 计算
    [ber(Idx),num(Idx),l(Idx)]=Receiver.Cal_BER(selectSignal);
    % 存储kk
    re_signal=[re_signal,selectSignal];
end


% CD compensation
paramEDC = struct();
paramEDC.L = param.Ltotal;
paramEDC.D = param.D;
paramEDC.Fc = 193.1e12;
paramEDC.Fs = fs;


Receiver.Nr.k=k;
% BER 计算
Receiver.Nr.nTrainSym=100;
% CD 补偿
paramEDC.D = param.D;
re_signal=cdc(re_signal,paramEDC);
re_signal=re_signal.';

[~,~,~,data_qam,~]=Receiver.Demodulation(re_signal);
% 转换为矩阵形式
data_mat=reshape(data_qam,nn.nModCarriers,[]);

% 重新提调制为ofdm
re_mod_signal=[];
for idx=1:k
    martix=data_mat(:,nn.nPkts*(idx-1)+1:idx*nn.nPkts);
    % 重新调制为ofdm
    ofdm_signal= nn.ofdm(martix);
    re_mod_signal=[re_mod_signal,ofdm_signal.'];
end
% 色散逆补偿
paramEDC.D = -param.D;
re_mod_signal=cdc(re_mod_signal,paramEDC);
re_mod_signal=re_mod_signal.';




alpha=1;
[ipd_error1]=iteraError((re_mod_signal),fs,alpha,Dc,Vdither(1),f1,f2);
[ipd_error2]=iteraError((signal),fs,alpha,Dc,Vdither(1),f1,f2);


figure;hold on;
plot(ipd_error1)
plot(ipd_error2)


% 误差归一化
Error_norm=pnorm(Error);
ipd_error_norm=pnorm(ipd_error1);

close all;
figure;hold on;
plot(ipd_error_norm)
plot(Error_norm)
legend('重建误差','理想误差')
title('归一化后的误差显示')


% 最小二乘法确定比例因子alpha

I=Error.*ipd_error1;
Q=ipd_error1.^2;

alpha_k=sum(I)/sum(Q)
