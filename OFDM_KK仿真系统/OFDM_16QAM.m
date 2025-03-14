clear;close all;
addpath('D:\001-处理中\相干算法\optical_communication')

OFDM_TX;

[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);
SpS=nn.Fs/nn.Rs;
%参考信号
ref_seq=reshape(qam_signal,1,[]);
ref_seq = repmat(ref_seq,1,100);
% 重复,达成数量限制
k=1;
% qam信号矩阵
ref_seq_mat=repmat(qam_signal,1,k);


signal=repmat(signal,k,1);
% dither 的频率
f1=40e3;
f2=60e3;
Fs_new=nn.Fs;
N=length(signal)/(Fs_new/f1);

%方案选择
type='ssb';


if strcmp(type,'dsb')
Vdither1 = Creat_dither(Fs_new,f1,N);
Vdither2 = Creat_dither(Fs_new,f2,N*(f2/f1));
elseif strcmp(type,'ssb')
% f1 and f2 的cos信号
Vdither1 = Creat_dither(Fs_new,f1,N);
Vdither_1=Creat_dither(Fs_new,f2,N*(f2/f1));
% f1 and f2 的sin信号
Vdither2 = Creat_dither1(Fs_new,f2,N*(f2/f1));
Vdither_2=Creat_dither1(Fs_new,f1,N);
end


%IQ
%RF and dither
% Vdither=0.3;
Amp=4;
% Dither 的幅度设置
V=0.1*9

%IQ
N=length(label);
% 相位偏差参数
W=240;
pilotIndex=1:1:W;
P_EST='none';
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

for Index=1:length(V)
    Vdither=V(Index);
%     snr=[5:1:25];
snr=10;
for index=1:length(snr)


Pi_dBm = 10;
Pi = 10^(Pi_dBm/10)*1e-3; %W

paramIQ.Vpi=9;
if strcmp(type,'dsb')
    paramIQ.VbI=-0.8*paramIQ.Vpi+Vdither*Vdither1;
    paramIQ.VbQ=-0.8*paramIQ.Vpi+Vdither*Vdither2;
elseif strcmp(type,'ssb')
    paramIQ.VbI=-0.8*paramIQ.Vpi+Vdither*(Vdither1+Vdither_1);
    paramIQ.VbQ=-0.8*paramIQ.Vpi+Vdither*(Vdither2+Vdither_2);
end
paramIQ.Vphi=paramIQ.Vpi/2;

Ai  = sqrt(Pi);
m=Modulation_index(Amp*signal.',paramIQ.Vpi,'ofdm');
fprintf('Modulation index=%3.3f\n',m);
sigTxo = iqm(Ai, Amp*signal, paramIQ);


%power
power=signalpower(sigTxo);
fprintf('optical signal power: %.2f dBm\n', 10 * log10(power / 1e-3));



paramPD.B =nn.Fs;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=nn.Fs;
%noise
% sigTxo=awgn(sigTxo,snr(index),'measured');
%pd
ipd_btb = pd(sigTxo, paramPD);

%% kk
nTrainSym =10*k;
rxsig = ipd_btb(1:2*floor(length(ipd_btb)/2)).';
c=0;
% c=V(index);
SSB_Sig_KK = KK_New(rxsig+c,nn.Fs,nn.Fs);


%下采样
data_kk = downsample(SSB_Sig_KK,nn.Nsam);  
% 解OFDM
data_kk_ofdm = reshape(data_kk,nn.fft_size+nn.nCP,[]);
data_kk_ofdm(1:nn.nCP,:) = [];
data_kk_ofdm = fft(data_kk_ofdm);
% get the modulated carriers
data_kk_ofdm = data_kk_ofdm(postiveCarrierIndex,:);
% channel estimation
rxTrainSymbol = data_kk_ofdm(:,1:nTrainSym);
qam_signal_mat=repmat(qam_signal,1,k);
refTrainSymbol = qam_signal_mat(:,1:nTrainSym);
Hf = mean(rxTrainSymbol./refTrainSymbol,2);
%     Hf = smooth(Hf,5);
% channel equalization
data_kk = data_kk_ofdm.*repmat(1./Hf,1,nn.nPkts*k);
% 相位偏差消除
if strcmp(P_EST,'phase')
phase_compensation;
end 

%保留信号矩阵
data_kk_mat=data_kk;
%归一化
 data_kk=data_kk(:);
 data_kk = data_kk./sqrt(mean(abs(data_kk(:)).^2));
% scatterplot(data_kk);
% figure;
% scatplot(real(data_kk),imag(data_kk),'circles');


     %% 解码，计算误码
    % 参考序列
    ref_seq_1 =qamdemod(ref_seq,nn.M,'OutputType','bit','UnitAveragePower',1);
    ref_seq_1=ref_seq_1(:);
    % 接收序列
    yyy = data_kk;

    yyy_1 = qamdemod(yyy,nn.M,'OutputType','bit','UnitAveragePower',1);
    yyy_1=yyy_1(:);

    [ber(Index,index),num(Index,index),error_location] = CalcBER(yyy_1,ref_seq_1); %计算误码率
    fprintf('Num of Errors = %d, BER = %1.7f\n',num(index),ber(index));

    Calc_BER_mat;
end
end



figure;
subplot(3,1,1)
plot(real(SSB_Sig_KK));hold on;
plot(imag(SSB_Sig_KK));
title('IQ输出')
subplot(3,1,2)
plot(ipd_btb)
title('PD接收')
subplot(3,1,3)
plot(1:idx,number_num);
title('误码数')
if 0
    figure;
semilogy(snr,ber)
xlabel('snr')
ylabel('ber')

end
%%
% % 生成路径
% dir_up='D:\PhD\Project\单边带光发射机自适应偏压控制\Simulation\传输_影响\BER_仿真';

% data_path = dir_up;
% if(~exist(data_path,'dir'))
%     mkdir(char(data_path));
% end
% save(data_path+"/BER_w_ssb_dither_10k_1%V.mat","ber");