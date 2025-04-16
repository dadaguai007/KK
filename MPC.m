% 复现：Minimum phase conditions in kramers-kronig optical receivers
% SSB 和 DSB dither 的 作用效果图

clear;clc;close all;
% addpath('D:\PhD\Codebase\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
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

flag_mon={'ssb','dsb'};
% flag_mon={'dsb'};
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
        leg_text=["Received signal";"Recovered signal using KK algorithm"];
        FontSize=14;
        figure;hold on;
        plot(real(s(1:2e5)),imag(s(1:2e5)),'ro',LineWidth=2);
        circles_plot(Ac,As,'MPC','Re','Im',[5 20],[-2*As 2*As],leg_text,FontSize)


        % 接收
        paramPD.B =f;
        paramPD.R =1;
        paramPD.type = 'ideal';
        paramPD.Fs=fs;
        % pd
        I_pd = pd(s, paramPD);
        
        % ReLU
        % 期望信号
%         xn=I_pd(1:2e6);
%         dn=As*exp(-1j*w*t+1j*pi/2);
%         sps=1;
%         taps_list=[15,9,0,0,0,0];
%         ref=8;
%         step_len=0.0001;
%         [y,e,w]=ReLU_dfe_lms(xn,dn,sps,ref,taps_list,step_len);

        % KK
        f_up=4*fs; % KK算法的采样率
        %     I_KK = KK(I_pd,fs,fs);
        [I_KK,ln_sig,H_sig] = KK_MPC(I_pd,fs,f_up);
        %下采样
        s_recovery=downsample(I_KK,f_up/fs);

    end

    % 绘图效果
    plot(real(s_recovery(1:4e5)),imag(s_recovery(1:4e5)),'bx',LineWidth=2);
    circles_plot(Ac,As,'MPC','Re','Im',[5 20],[-2*As 2*As],leg_text,FontSize)

    if strcmp(type,'ssb')
        % 变量名
        [IfdBm_Dith_ssb , ~ ]=mon_ESA_flag(ln_sig,f_up,0);
        [ IfdBm_Dith_ssb_Ht, Fre ]=mon_ESA_flag(H_sig,f_up,0);
        [ IfdBm_Dith_ssb_receiver, fre ]=mon_ESA_flag(s_recovery,fs,0);
         [ IfdBm_Dith_ssb_I, fre ]=mon_ESA_flag(I_pd,fs,0);
    elseif strcmp(type,'dsb')
        [IfdBm_Dith_dsb , ~ ]=mon_ESA_flag(ln_sig,f_up,0);
        [ IfdBm_Dith_dsb_Ht, Fre ]=mon_ESA_flag(H_sig,f_up,0);
         [ IfdBm_Dith_dsb_receiver, fre ]=mon_ESA_flag(s_recovery,fs,0);
         [ IfdBm_Dith_dsb_I, fre ]=mon_ESA_flag(I_pd,fs,0);
    end
   WB.updata(index);
end
WB.closeWaitBar();


figure;hold on;
plot(Fre,IfdBm_Dith_ssb,'r');
plot(Fre,IfdBm_Dith_dsb,'b--');
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('ln_{signal}')
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 14);
set(gcf,'Position', [0, 0, 480, 400]);
set(gca, 'LineWidth', 1.25);
set(gca,'XLim',[-f_up/2/1e9 f_up/2/1e9],'YLim',[-50 30]);



figure;hold on;
plot(Fre,IfdBm_Dith_ssb_Ht,'r');
plot(Fre,IfdBm_Dith_dsb_Ht,'b--');
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('HT_{signal}')
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 14);
set(gcf,'Position', [0, 0, 480, 400]);
set(gca, 'LineWidth', 1.25);
set(gca,'XLim',[-f_up/2/1e9 f_up/2/1e9],'YLim',[-50 30]);


color=distinguishable_colors(20);
figure;hold on;
plot(fre,IfdBm_Dith_dsb_receiver,'Color',color(2,:));
plot(fre,IfdBm_Dith_ssb_receiver,'Color',color(12,:));
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('After KK')
legend('DSB','SSB');
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 14);
set(gcf,'Position', [0, 0, 480, 400]);
set(gca, 'LineWidth', 1.25);
set(gca,'XLim',[-fs/2/1e9 fs/2/1e9],'YLim',[-80 30]);


figure;
plot(fre,IfdBm_Dith_dsb_receiver,'Color',color(12,:));
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('After KK')
legend('DSB');
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 14);
set(gcf,'Position', [0, 0, 480, 400]);
set(gca, 'LineWidth', 1.25);
set(gca,'XLim',[-fs/2/1e9 fs/2/1e9],'YLim',[-400 50]);


figure;
plot(fre,IfdBm_Dith_ssb_receiver,'Color',color(12,:));
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('After KK')
legend('SSB');
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 14);
set(gcf,'Position', [0, 0, 480, 400]);
set(gca, 'LineWidth', 1.25);
set(gca,'XLim',[-fs/2/1e9 fs/2/1e9],'YLim',[-400 50]);


figure;
plot(fre,IfdBm_Dith_dsb_I,'Color',color(12,:));
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('After PD')
legend('DSB');
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 14);
set(gcf,'Position', [0, 0, 480, 400]);
set(gca, 'LineWidth', 1.25);
set(gca,'XLim',[-fs/2/1e9 fs/2/1e9],'YLim',[-300 50]);


figure;
plot(fre,IfdBm_Dith_ssb_I,'Color',color(12,:));
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('After PD')
legend('SSB');
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 14);
set(gcf,'Position', [0, 0, 480, 400]);
set(gca, 'LineWidth', 1.25);
set(gca,'XLim',[-fs/2/1e9 fs/2/1e9],'YLim',[-300 50]);

