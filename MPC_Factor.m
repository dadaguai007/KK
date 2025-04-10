% 复现：DSB MPC 影响因素 探究

clear;clc;close all;
addpath('D:\PhD\Codebase\')
% addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
addpath('Fncs\')

% 载波幅度
Ac=4;
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
        % 无扰动信号
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



        % 接收
        paramPD.B =f;
        paramPD.R =1;
        paramPD.type = 'ideal';
        paramPD.Fs=fs;
        % pd
        I_pd = pd(s, paramPD);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        %  无扰动信号/ 试验扰动影响       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 无扰动项
        I_pd_Tone = pd(SS.', paramPD);
        % 平方
        S_dither=dither1+1j*dither2;
        I_pd_dither = pd(S_dither.', paramPD);

        % 平方项
        I_pd_Total=I_pd_Tone+I_pd_dither;  %加上平方项

        % 信号交调项

        S=As*exp(-1j*w*t+1j*pi/2).*dither1+conj(As*exp(-1j*w*t+1j*pi/2)).*dither1;

        S1=real(As*exp(-1j*w*t+1j*pi/2).*conj(1j*dither2)+conj(As*exp(-1j*w*t+1j*pi/2)).*(1j*dither2));
        % 平方项+拍频
        I_pd_Total1=I_pd_Total+S.';

        % 平方项+ 两个拍频项
        I_pd_Total2=I_pd_Total1+S1.';

        % dither的扩大
        S2=S_dither*Ac+Ac*conj(S_dither);

        % dither扩大项 + 两个拍频项
        I_pd_Total3=I_pd_Tone+S2.'+S.'+S1.';

        % dither扩大项 + 一个拍频项
        I_pd_Total4=I_pd_Tone+S2.'+S.';
        % dither扩大项
        I_pd_Total5=I_pd_Tone+S2.';
        % 拍频项
        I_pd_Total6=I_pd_Tone+S1.'+S.';
        % 平方项
        I_pd_Total7=I_pd_Total+S2.';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        % KK算法            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        f_up=4*fs; % KK算法的采样率
        
        [I_KK,~,~] = KK_MPC(I_pd_Total3,fs,f_up);
        [I_KK1,~,~] = KK_MPC(I_pd_Total4,fs,f_up);
        [I_KK2,~,~] = KK_MPC(I_pd_Total5,fs,f_up);
        [I_KK3,~,~] = KK_MPC(I_pd_Total6,fs,f_up);
        [I_KK4,~,~] = KK_MPC(I_pd_Total,fs,f_up);
        [I_KK5,~,~] = KK_MPC(I_pd_Total2,fs,f_up);
        [I_KK6,~,~] = KK_MPC(I_pd_Total7,fs,f_up);
        [I_KK_Total,~,~] = KK_MPC(I_pd,fs,f_up);
        
        %下采样
        s_recovery=downsample(I_KK,f_up/fs);
        s_recovery1=downsample(I_KK1,f_up/fs);
        s_recovery2=downsample(I_KK2,f_up/fs);
        s_recovery3=downsample(I_KK3,f_up/fs);
        s_recovery4=downsample(I_KK4,f_up/fs);
        s_recovery5=downsample(I_KK5,f_up/fs);
        s_recovery6=downsample(I_KK6,f_up/fs);
        s_recovery_total=downsample(I_KK_Total,f_up/fs);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      绘图效果            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FontSize=14;
    leg_text=["Received signal";"Recovered signal using KK algorithm"];
    figure;hold on;
    plot(real(s(1:2e5)),imag(s(1:2e5)),'ro',LineWidth=2);
    plot(real(s_recovery(1:2e5)),imag(s_recovery(1:2e5)),'bx',LineWidth=2);
    circles_plot(Ac,As,'载波-扰动拍频+两个扰动-扰动拍频','Re','Im',[5 20],[-2*As 2*As],leg_text,FontSize)
    

    figure;hold on;
    plot(real(s(1:2e5)),imag(s(1:2e5)),'ro',LineWidth=2);
    plot(real(s_recovery1(1:2e5)),imag(s_recovery1(1:2e5)),'bx',LineWidth=2);
    circles_plot(Ac,As,'载波-扰动拍频+一个扰动-扰动拍频','Re','Im',[5 20],[-2*As 2*As],leg_text,FontSize)

    figure;hold on;
    plot(real(s(1:2e5)),imag(s(1:2e5)),'ro',LineWidth=2);
    plot(real(s_recovery2(1:2e5)),imag(s_recovery2(1:2e5)),'bx',LineWidth=2);
    circles_plot(Ac,As,'载波-扰动拍频','Re','Im',[5 20],[-2*As 2*As],leg_text,FontSize)

    figure;hold on;
    plot(real(s(1:2e5)),imag(s(1:2e5)),'ro',LineWidth=2);
    plot(real(s_recovery3(1:2e5)),imag(s_recovery3(1:2e5)),'bx',LineWidth=2);
    circles_plot(Ac,As,'两个扰动-扰动拍频','Re','Im',[5 20],[-2*As 2*As],leg_text,FontSize)

     figure;hold on;
    plot(real(s(1:2e5)),imag(s(1:2e5)),'ro',LineWidth=2);
    plot(real(s_recovery4(1:2e5)),imag(s_recovery4(1:2e5)),'bx',LineWidth=2);
    circles_plot(Ac,As,'扰动平方','Re','Im',[5 20],[-2*As 2*As],leg_text,FontSize)
    
    figure;hold on;
    plot(real(s(1:2e5)),imag(s(1:2e5)),'ro',LineWidth=2);
    plot(real(s_recovery5(1:2e5)),imag(s_recovery5(1:2e5)),'bx',LineWidth=2);
    circles_plot(Ac,As,'扰动平方+两个扰动-扰动拍频','Re','Im',[5 20],[-2*As 2*As],leg_text,FontSize)

    figure;hold on;
    plot(real(s(1:2e5)),imag(s(1:2e5)),'ro',LineWidth=2);
    plot(real(s_recovery5(1:2e5)),imag(s_recovery5(1:2e5)),'bx',LineWidth=2);
    circles_plot(Ac,As,'扰动平方+载波-扰动拍频','Re','Im',[5 20],[-2*As 2*As],leg_text,FontSize)


    figure;hold on;
    plot(real(s(1:2e5)),imag(s(1:2e5)),'ro',LineWidth=2);
    plot(real(s_recovery_total(1:2e5)),imag(s_recovery_total(1:2e5)),'bx',LineWidth=2);
    circles_plot(Ac,As,'所有扰动项','Re','Im',[5 20],[-2*As 2*As],leg_text,FontSize)

    WB.updata(index);
end
WB.closeWaitBar();


