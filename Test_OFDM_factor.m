% dither 对 多载波信号的影响
clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Project\Base_Code\Base\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
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
k=30;
% qam信号矩阵
ref_seq_mat=repmat(qam_signal,1,k);

% 信号复制
signal=repmat(signal,k,1);
% dither 的频率处理
f1=400e6;
f2=600e6;
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

VbI_sin = Vdither(1)*Creat_dither1(Fs_new,f1,N);
VbQ_sin = Vdither(1)*Creat_dither1(Fs_new,f2,N*(f2/f1));
% 验证 信号是否正确
a1=Vdither(1)/2*Creat_ssb(Fs_new,f1,N);
a2=Vdither(1)/2*Creat_ssb1(Fs_new,f1,N);

b1=1j*Vdither(1)/2*Creat_ssb(Fs_new,f2,N*(f2/f1));
b2=1j*Vdither(1)/2*Creat_ssb1(Fs_new,f2,N*(f2/f1));

error_tone=VbI-(a1+a2);
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
% 信号拍频项
E=pd(Ac+signal,paramPD);
% dither 平方项
ipd_pingfang = pd(Dither, paramPD);
% dither的扩大
S2=Dither*Ac+Ac*conj(Dither);
S2_1=Ac*VbI+Ac*conj(VbI);
S2_2=Ac*1j*VbQ+Ac*conj(1j*VbQ);
% 拍频项 1
S=As*signal.*VbI+conj(As*signal).*VbI;
% 拍频项 2
S1=real(As*signal.*conj(1j*VbQ)+conj(As*signal).*(1j*VbQ));

% 误差，探寻分量设置是否正确
error=ipd_btb-ipd_pingfang-S2-S-S1-E;


% I 路
% 负频率 dither

% 与载波拍频
E2_S2=Ac*VbI;
% 与信号拍频
E2_S=signal.*conj(a1)+conj(signal).*a1;
E2_S_hat=real(signal).*VbI-imag(signal).*VbI_sin;

% 验证 负dither 与 信号拍频的简便表达形式
error_ff=E2_S-E2_S_hat;

% 自拍频
E2_self=a1.*conj(a1);

% 正频率 dither

% 与载波拍频
E1_S2=Ac*VbI;
% 与信号拍频
E1_S=signal.*conj(a2)+conj(signal).*a2;
E1_S_hat=real(signal).*VbI+imag(signal).*VbI_sin;

% 验证 正dither 与 信号拍频的简便表达形式
error_ff1=E1_S-E1_S_hat;
% 自拍频
E1_self=a2.*conj(a2);

% 验证 dither 信号拍频是否正确
error_S=S-E2_S-E1_S;
% 验证 dither 载波拍频是否正确
error_S2=S2_1-E2_S2-E1_S2;


% Q 路
% 负频率

% 与载波拍频
E2_S1=Ac*VbQ_sin;
% 与信号拍频
E2_S_car=signal.*conj(b1)+conj(signal).*b1;
E2_S_car_hat=real(signal).*VbQ_sin+imag(signal).*VbQ;


% 正频率

% 与载波拍频
E1_S1=-Ac*VbQ_sin;
% 与信号拍频
E1_S_car=signal.*conj(b2)+conj(signal).*b2;
E1_S_car_hat=-real(signal).*VbQ_sin+imag(signal).*VbQ;

% 验证 Q 路 dither 信号拍频是否正确
error_S1=S1-E2_S_car_hat-E1_S_car_hat;


% 接收信号E2，E1，E
% 信号
E=2*real(signal)*Ac+Ac.^2;
% 信号自拍频
E_ssbi=signal.*conj(signal);
% 正频率
E1=E1_S_hat+E1_S_car_hat;
% 负频率
E2=E2_S_hat+E2_S_car_hat;

% 定义pwelch函数的参数
window = hamming(256); % 使用汉明窗
noverlap = 128; % 重叠128个样本
nfft = 1024; % FFT点数

% 计算每个信号的功率谱密度
[PxxE,F1] = pwelch(E,window,noverlap,nfft,fs);
[Pxxssbi,F2] = pwelch(E_ssbi,window,noverlap,nfft,fs);
[PxxE1,F3] = pwelch(E1,window,noverlap,nfft,fs);
[PxxE2,F4] = pwelch(E2,window,noverlap,nfft,fs);
[Pxxpd,F5] = pwelch(ipd_btb,window,noverlap,nfft,fs);


figure;
plot(F1,10*log10(PxxE),'LineWidth',2);
hold on;
plot(F2,10*log10(Pxxssbi),'LineWidth',2);
plot(F3,10*log10(PxxE1),'LineWidth',2);
plot(F4,10*log10(PxxE2),'LineWidth',2);
plot(F5,10*log10(Pxxpd),'LineWidth',2);
legend('信号','ssbi','E1-正频率','E2-负频率','PD接收')