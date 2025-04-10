clear;clc;close all;
% addpath('D:\PhD\Codebase\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
addpath('Fncs\')
% 信号生成
Muti_Tone_Tx;
% dither
f1=40e3;
f2=60e3;
w_dither1=2*pi*f1;
w_dither2=2*pi*f2;


% 载波幅度
Ac=8;


flag_mon={'ssb','dsb'};
for i=1:length(flag_mon)
% 信号幅的大小
Ratio_power=0.2;

% dither 信号 选择不同的dither方案
type=char(flag_mon(i));

ratio_power=Ratio_power;
if strcmp(type,'dsb')
    dither1=ratio_power*As*cos(w_dither1*t_up);
    dither2=1j*ratio_power*As*cos(w_dither2*t_up);
elseif strcmp(type,'ssb')
    dither1=ratio_power*As*exp(1j*w_dither1*t_up);
    dither2=ratio_power*As*exp(1j*w_dither2*t_up);
end
% 传输信号
Sig_TX=signal+dither1+dither2+Ac;

% 纯净信号

Sig_=signal+Ac;

% 除去载波
signal_power=signalpower(Sig_TX-Ac);
% Cal CSPR
Carry_power=Ac.^2;
CSPR = Carry_power/signal_power;
CSPR=10*log10(CSPR);
CSPR=round(CSPR,2);
fprintf("the CSPR is %1.2fdB\n",CSPR);

% 接收
paramPD.B =fs;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=fs;
% pd
I_pd = pd(Sig_TX, paramPD);
btb=pd(Sig_,paramPD);

% 纯净信号
Signal_Beat(i,:)=btb;
% 频谱存储
[ IfdBm, Fre]=mon_ESA_flag((Signal_Beat(i,:)),fs,0);
Signal_Beat_PowerdBm(i,:)=IfdBm;

% 接收信号
Signal(i,:)=I_pd;
% 频谱存储
[ IfdBm, ~ ]=mon_ESA_flag((Signal(i,:)),fs,0);
Signal_PowerdBm(i,:)=IfdBm;

end
color=distinguishable_colors(20);
% 单边带
figure;hold on;
plot(Fre,Signal_PowerdBm(1,:),'Color',color(12,:));
plot(Fre,Signal_Beat_PowerdBm(1,:),'Color',color(2,:));

FontSize=14;
flag=struct();
flag.LegendON_OFF=1;
flag.Legendflage=0;
% 标题
LegendArray=["Noise Signal","Signal-Beat"];

Plotter('Receiver Signal Spectrum','Frequency (GHz)','Magnitude (dB)',[-fs/2/1e9 fs/2/1e9],...
    [-250 20],LegendArray,flag,FontSize);


% 双边带
figure;hold on;
plot(Fre,Signal_PowerdBm(2,:),'Color',color(12,:));
plot(Fre,Signal_Beat_PowerdBm(2,:),'Color',color(2,:));

FontSize=14;
flag=struct();
flag.LegendON_OFF=1;
flag.Legendflage=0;
% 标题
LegendArray=["Noise Signal","Signal-Beat"];

Plotter('Receiver Signal Spectrum','Frequency (GHz)','Magnitude (dB)',[-fs/2/1e9 fs/2/1e9],...
    [-250 20],LegendArray,flag,FontSize);