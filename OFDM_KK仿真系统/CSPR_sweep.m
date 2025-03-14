clear;close all;clc;
addpath('Fncs\')
% addpath('D:\001-处理中\相干算法\optical_communication')
OFDM_TX;
% 生成信号
[y1,y2,signal,qam_signal,postiveCarrierIndex]=nn.Output();
label = nn.ofdm(qam_signal);
SpS=nn.Fs/nn.Rs;

scale_factor = max(max(abs(real(signal))),max(abs(imag(signal))));
signal = signal./scale_factor;

%参考信号
ref_seq=reshape(qam_signal,1,[]);
ref_seq = repmat(ref_seq,1,100);
% 重复信号
k=1;
% qam信号矩阵
ref_seq_mat=repmat(qam_signal,1,k);


signal=repmat(signal,k,1);
% dither 的频率处理
f1=40e3;
f2=60e3;
Fs_new=nn.Fs;
N=length(signal)/(Fs_new/f1);

%方案选择

flagtag=0;% 绘图标致


m1=0.1135;

%IQ
%RF and dither
Amp=1.7;%1.4 1.7 2 2.75 3.4 4.1
% Dither 的幅度设置
alfa=0;
V=alfa*9;

%IQ
N=length(label);
% 相位均衡参数
W=nn.nModCarriers;
pilotIndex=1:1:W;
P_EST='none';
if 0
    CR =2.3; % Clipping Ratio
    signal_orgin=signal;
    signal = clipping(signal.',CR);
else
    signal=signal.';
end

for Index=1:length(V)
    Vdither=V(Index);
    dx=0.01;
    beta=0:dx:1;
    for index=1:length(beta)

        Pi_dBm = 10;
        Pi = 10^(Pi_dBm/10)*1e-3; %W

        paramIQ.Vpi=9;
        beta1=beta(index);

        paramIQ.VbI=beta1*paramIQ.Vpi;
        paramIQ.VbQ=beta1*paramIQ.Vpi;

        paramIQ.Vphi=paramIQ.Vpi/2;

        Ai  = sqrt(Pi);
        m=Modulation_index(Amp*signal.',paramIQ.Vpi,'ofdm');
        fprintf('Modulation index=%3.3f\n',m);
        sigTxo = iqm(Ai, Amp*signal, paramIQ);


        %power
        power=signalpower(sigTxo);
        fprintf('optical signal power: %.2f dBm\n', 10 * log10(power / 1e-3));

        sigRxo=sigTxo;

        %PD
        paramPD.B =nn.Fs;
        paramPD.R =1;
        paramPD.type = 'ideal';
        paramPD.Fs=nn.Fs;
        %pd
        ipd_btb = pd(sigRxo, paramPD);

        % CSPR theory
        paramIQ.Vpi=9;
        paramIQ.VbI=beta1*paramIQ.Vpi;
        paramIQ.VbQ=beta1*paramIQ.Vpi;
        paramIQ.Vphi=paramIQ.Vpi/2;
        if 1
            Eout1 = iqm(Ai, 0, paramIQ);
            Eout_s=sigTxo-Eout1*ones(1,length(sigTxo));
            power1=signalpower(Eout1);
            Ps=signalpower(Eout_s);
            CSPR = power1/(Ps);
            CSPR=real(10*log10(CSPR));
        else
            U=mean(ipd_btb);
            V=var(ipd_btb);
            x1=sqrt(U.^2-V);
            x2=U-x1;
            CSPR = x1/x2;
            CSPR=10*log10(CSPR);

        end
        fprintf('CSPR1 = %1.7f\n',CSPR);
        CSPR_theory(index)=(CSPR);

        % CSPR sim
        Vpi=paramIQ.Vpi;
        c=(pi.^2)/(4*(Vpi.^2));
        d=1/2*(1-cos(beta1*pi));
        d1=1/2*(1+cos(beta1*pi));
        Mindx = str2double(sprintf('%.3f', m/2));
        fprintf('m = %1.2f\n',Mindx);
        %         p1=((m1*Vpi).^2)/9; 自己设置的调制深度
        p2=((Mindx*Vpi).^2)/9;
        p=var(Amp*real(signal));
        CSPR4 = d1/(p2*c*d) ;
        CSPR4=10*log10(CSPR4);

        fprintf('CSPR4 is %1.3f\n',CSPR4);
        CSPR_sim(index)=CSPR4;
        
        % CSPR theory
        Eout_10 = MZM_SingleDrive(ones(1,length(signal)),Amp*real(signal),paramIQ.VbI,Vpi,1);
        Eout_9 = MZM_SingleDrive(ones(1,length(signal)),zeros(1,length(signal)),paramIQ.VbI,Vpi,1);
        E_signal=Eout_10-Eout_9;
        power_total_signal=signalpower(E_signal);
        power_signal=signalpower(Amp*real(signal));

        %CSPR5=d1/(power_signal*c*d);
        CSPR5=d1/(power_total_signal);
        CSPR5=10*log10(CSPR5);
        CSPR_theory1(index)=real(CSPR5);

        POWER(index)=power_total_signal;
        POWER_SIG(index)=power_signal*c*d;
    end
end

figure;hold on;
plot((POWER))
plot(POWER_SIG)
title('信号功率对比')
figure;hold on;
plot(CSPR_theory1);
plot(CSPR_sim)
% plot(CSPR_theory);
title('CSPR')

figure;hold on;
plot(CSPR_theory1(1:end-1));
plot(CSPR_theory(1:end-1));
title('MZM vs IQ CSPR')

if 0
    Title='CSPR_theory1';
    datapath = strcat('Output\',Title);
    if ~exist(datapath,'dir')
        mkdir(datapath);
    end
    save(sprintf('%s\\theory_%.2f.mat',datapath,Mindx),'CSPR_theory1');


    Title='CSPR_sim';
    datapath1 = strcat('Output\',Title);
    if ~exist(datapath1,'dir')
        mkdir(datapath1);
    end
    save(sprintf('%s\\sim_%.2f.mat',datapath1,Mindx),'CSPR_sim');
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