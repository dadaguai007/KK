% dither 对 多载波信号的影响
clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Project\Base_Code\Base\')
% addpath('D:\BIT_PhD\Base_Code\Codebase_using\')
addpath("Plot\")

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

%  全部解码
if 1
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
    Vpi=9;
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
        10*k, ...            % 信道训练长度
        1:1:ofdmPHY.nModCarriers, ...    %导频位置
        1, ...             % 选取第一段信号
        ref_seq, ...       % 参考序列
        ref_seq_mat, ...    % qam 矩阵
        'off', ...         % 是否采用CPE
        'off', ...         % 对所有载波进行相位补偿
        'KK');             % 接收方式

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

        % 信号预处理
        [ReceivedSignal,dc]=Receiver.Total_Preprocessed_signal(pd_receiver);
        Dc=mean(ReceivedSignal);
        % BER 计算
        [ber_total(index),num_total(index)]=Receiver.Cal_BER(ReceivedSignal);
        if 1
            disp('进入循环');
            % 设置单边带信号
            VbQ_sin = Vdither(1)*Creat_dither1(Fs_new,f2,N*(f2/f1));
            VbI_sin=Vdither(1)*Creat_dither1(Fs_new,f1,N);
            
            % 赋值信号进入迭代
            ipd_pd=pd_receiver;
            re=ReceivedSignal;
            for jj=1:20

                dd=real(Dc);
                % 重新计算E2和E1
                Rece_remove_dc=re-dd;

                % 载波与dither 拍频
                E5=real(dd)*VbI+real(dd)*VbQ_sin;
                E4=real(dd)*VbI-real(dd)*VbQ_sin;
                E1=E5+E4;
                % 信号与dither拍频
                %         I_beat=Rece_remove_dc.*conj(VbI_ssb)+conj(Rece_remove_dc).*VbI_ssb;
                %         Q_beat=Rece_remove_dc.*conj(VbQ_ssb)+conj(Rece_remove_dc).*VbQ_ssb;

                %  负频率
                I_beat=real(Rece_remove_dc).*VbI-imag(Rece_remove_dc).*VbI_sin;
                Q_beat=real(Rece_remove_dc).*VbQ_sin+imag(Rece_remove_dc).*VbQ;

                E2=I_beat+Q_beat;
                % 正频率
                %         I_beat1=real(Rece_remove_dc).*VbI+imag(Rece_remove_dc).*VbI_sin;
                %         Q_beat1=-real(Rece_remove_dc).*VbQ_sin+imag(Rece_remove_dc).*VbQ;
                %
                %         E3=I_beat1+Q_beat1;


                alpha=0.02;
                ipd_error=alpha*(E1+E2);

                ipd_pd=ipd_pd-ipd_error;




                if 0
                    re_signal=[];
                    % 发射机参数
                    ofdmPHY=nn;
                    for i=1:k

                        %%---------------------------------------        解码       ---------------------------%%
                        Receiver=OFDMreceiver( ...
                            ofdmPHY, ...       %%% 发射机传输的参数
                            ofdmPHY.Fs, ...    %   采样
                            6*ofdmPHY.Fs, ...  % 上采样
                            ofdmPHY.nPkts, ...            % 信道训练长度
                            1:1:ofdmPHY.nModCarriers, ...    %导频位置
                            i, ...             % 选取第一段信号
                            ref_seq, ...       % 参考序列
                            qam_signal, ...    % qam 矩阵
                            'off', ...         % 是否采用CPE
                            'off', ...         % 对所有载波进行相位补偿
                            'KK');             % 接收方式

                        % 信号预处理
                        [receive,~]=Receiver.Preprocessed_signal(ipd_pd);
                        % BER 计算
                        [ber,num(i),l(i)]=Receiver.Cal_BER(receive);
                        re_signal=[re_signal,receive];

                    end

                    fprintf('分组解码的BER = %1.7f\n',sum(num)/sum(l));

                    % 用分组解码出来的信号代替
                    re=re_signal;
                    Dc=mean(re_signal);
                    Receiver.Total_Preprocessed_signal(ipd_pd);
                end
%                 ber_group(jj)=sum(num)/sum(l);


                % 信号预处理
                [re,~]=Receiver.Total_Preprocessed_signal(ipd_pd);
                Dc=mean(re);
                Receiver.Nr.nTrainSym = ofdmPHY.nPkts;
                % BER 计算
                [ber_total1(jj),num_total1]=Receiver.Cal_BER(re);
            end
            ber_total_iter(index)=min(ber_total1);
%             ber_total_iter_Group(index)=min(ber_group);
        end
            WB.updata(index);
        
    end
    WB.closeWaitBar();
end





berplot = BERPlot_David();
% 间隔
berplot.interval=3;
% 字号
berplot.Config.FontSize = 14;
berplot.flagThreshold=1;
berplot.flagRedraw=0;
berplot.flagAddLegend=1;
% BER=[ber_total.';ber_total_iter;ber_total_iter_Group];
BER=[ber_total.';ber_total_iter;];
LengendArrary=["80km w/o ","80km w SIC"];
berplot.multiplot(Eb_N0_dB,BER,LengendArrary);