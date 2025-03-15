% 复现：ReLU 效果 (无明显效果)
clear;clc;close all;
% addpath('D:\PhD\Codebase\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')

addpath('Fncs\')
% 载波幅度
Ac=5;
% 信号幅度
As=1;
% 频率
f=1e9;
f1=40e3;
f2=60e3;
SpS=4;
% 采样率
fs=f*SpS;
% Time array
T=1/fs;
t=0:T:1e-3-T;
% 角频率
w=2*pi*f;
w_dither1=2*pi*f1;
w_dither2=2*pi*f2;


flag_mon={'dsb'};
WB = OCG_WaitBar(length(flag_mon));
for index= 1:length(flag_mon)

    % dither 信号 选择不同的dither方案
    type=char(flag_mon(index));
    % 信号幅的大小
    Ratio_power=0.1;
    % Ratio_power=0:0.01:0.1;
    for i=1:length(Ratio_power)
        ratio_power=Ratio_power(i);
        if strcmp(type,'dsb')
            dither1=ratio_power*As*cos(w_dither1*t);
            dither2=ratio_power*As*cos(w_dither2*t);
        elseif strcmp(type,'ssb')
            dither1=ratio_power*As*exp(-1j*w_dither1*t);
            dither2=ratio_power*As*exp(-1j*w_dither2*t);
        end
        % 信号
        if strcmp(type,'dsb')
            s=Ac+As*exp(-1j*w*t+1j*pi/2)+dither1+1j*dither2;
            s=s.';
            signal=As*exp(-1j*w*t)+dither1+1j*dither2;
        else
            s=Ac+As*exp(-1j*w*t+1j*pi/2)+dither1+dither2;
            s=s.';
            signal=As*exp(-1j*w*t)+dither1+dither2;
        end
        SS=Ac+As*exp(-1j*w*t+1j*pi/2);
        signal_power=signalpower(signal);

        % Cal CSPR
        Carry_power=Ac.^2;
        CSPR = Carry_power/signal_power;
        CSPR=10*log10(CSPR);
        CSPR=round(CSPR,2);
        text='CSPR=';
        leg_text=strcat(text,num2str(CSPR),'dB');
        fprintf("the CSPR is %1.2fdB\n",CSPR);

        % 绘图效果
        FontSize=14;
        figure;hold on;
        plot(real(s(1:2e5)),imag(s(1:2e5)),'ro',LineWidth=2);
        circles_plot(Ac,As,'MPC','Real','Imag',[5 20],[-2*As 2*As],leg_text,FontSize)


        % 接收
        paramPD.B =f;
        paramPD.R =1;
        paramPD.type = 'ideal';
        paramPD.Fs=fs;
        % pd
        I_pd = pd(s, paramPD);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        % ReLU算法         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %         期望信号
        xn=I_pd(1:1e6);
        d=Ac+As*exp(-1j*w*t+1j*pi/2);
        dn=pd(d,paramPD);
        sps=1;
        taps_list=[11,3,0,0,0,0];
        ref=1;
        step_len=0.0001;
        [y,e,mu]=ReLU_dfe_lms(xn,dn,sps,ref,taps_list,step_len);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        % KK算法            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        f_up=4*fs; % KK算法的采样率
        %     I_KK = KK(I_pd,fs,fs);
        [I_KK,ln_sig,H_sig] = KK_MPC(I_pd,fs,f_up);
        %下采样
        s_recovery=downsample(I_KK,f_up/fs);

    end

    % 绘图效果
    plot(real(s_recovery(1:2e5)),imag(s_recovery(1:2e5)),'bx',LineWidth=2);
    circles_plot(Ac,As,'MPC','Real','Imag',[5 20],[-2*As 2*As],leg_text,FontSize)

    WB.updata(index);
end
WB.closeWaitBar();


mon_ESA(y,fs);
mon_ESA(I_pd,fs);