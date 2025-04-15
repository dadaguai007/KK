% dither 对 多载波信号的影响
clear;close all;clc;
addpath('Fncs\')
addpath('D:\PhD\Project\Base_Code\Base\')
% addpath('D:\BIT_PhD\Base_Code\Codebase_using\')

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


%IQ
N=length(label);

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

    phi=0.87; % Vbias 偏移程度
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
    Eb_N0_dB=10:30;
    ber_total=zeros(length(Eb_N0_dB),1);
    num_total=zeros(length(Eb_N0_dB),1);
    WB = OCG_WaitBar(length(Eb_N0_dB));
    for index=1:length(Eb_N0_dB)
        % 噪声
        noise=EbN0_dB(sigRxo,Eb_N0_dB(index));
        % 加入噪声
        pd_receiver = pnorm(ipd_btb)+(noise);
        % 信号预处理
        [ReceivedSignal,~]=Receiver.Total_Preprocessed_signal(pd_receiver);
        % 测试两个DC的取值
        Dc=mean(ReceivedSignal);
        % BER 计算
        [ber_total(index),num_total(index)]=Receiver.Cal_BER(ReceivedSignal);

        Error=[];
        for i=1:k
            % 选取某段进行解码(可进行优化)
            %%---------------------------------------        解码       ---------------------------%%
            Receiver=OFDMreceiver( ...
                ofdmPHY, ...       %%% 发射机传输的参数
                ofdmPHY.Fs, ...    %   采样
                6*ofdmPHY.Fs, ...  % 上采样
                ofdmPHY.nPkts, ...            % 信道训练长度
                1:1:ofdmPHY.nModCarriers, ...    %导频位置
                i, ...             % 选取一段信号
                ref_seq, ...       % 参考序列
                qam_signal, ...    % qam 矩阵
                'off', ...         % 是否采用CPE
                'off', ...         % 对所有载波进行相位补偿
                'KK');             % 接收方式
            % obj.Nr.squ_num 选取第n段。
%             Receiver.Button.Clipping='on';
%             Receiver.Nr.CL=0.2;
            % 信号预处理
            [receive,dc]=Receiver.Preprocessed_signal(pd_receiver);
            [signal_ofdm_martix,data_ofdm_martix,Hf,data_qam,qam_bit]=Receiver.Demodulation(receive);

            % BER 计算
            [ber,num]=Receiver.Cal_BER(receive);

            % 选取性能较好段，进行重新调制
            [ofdm_signal,~] = nn.ofdm(data_ofdm_martix);
            % 补上直流
            Re_Signal=Dc+ofdm_signal;

            rr=receive.'-Re_Signal;
            % 计算误差
            Error=cat(1,Error,rr);

        end

        Signal=ReceivedSignal-Error.';
        % 更正解码参数
        Receiver.Total_Preprocessed_signal(ipd_btb);
        %解码
        [ber_total1(index),num_total1(index)]=Receiver.Cal_BER(Signal);


        % 再进行一次循环
        Error1=[];
        for j=1:k
            Receiver=OFDMreceiver( ...
                ofdmPHY, ...       %%% 发射机传输的参数
                ofdmPHY.Fs, ...    %   采样
                6*ofdmPHY.Fs, ...  % 上采样
                ofdmPHY.nPkts, ...            % 信道训练长度
                1:1:ofdmPHY.nModCarriers, ...    %导频位置
                j, ...             % 选取一段信号
                ref_seq, ...       % 参考序列
                qam_signal, ...    % qam 矩阵
                'off', ...         % 是否采用CPE
                'off', ...         % 对所有载波进行相位补偿
                'KK');             % 接收方式
            Y=Signal(Receiver.ofdmPHY.len*(Receiver.Nr.squ_num-1)+1:Receiver.ofdmPHY.len*Receiver.Nr.squ_num);
            [signal_ofdm_martix_RE,data_ofdm_martix_RE,Hf_RE,data_qam_RE,qam_bit_RE]=Receiver.Demodulation(Y);
            [ber,num]=Receiver.Direcct_Cal_BER(data_qam_RE(:));

            % 选取性能较好段，进行重新调制
            [ofdm_signal1,~] = nn.ofdm(data_ofdm_martix_RE);
            % 补上直流
            Re_Signal1=Dc+ofdm_signal1;

            rr1=Y.'-Re_Signal1;
            % 计算误差
            Error1=cat(1,Error1,rr1);

        end
        Signal1=Signal-Error1.';
        % 更正解码参数
        Receiver.Total_Preprocessed_signal(ipd_btb);
        %解码
        [ber_total2(index),num_total2(index)]=Receiver.Cal_BER(Signal1);



        % 再进行一次循环
        Error2=[];
        for g=1:k
            Receiver=OFDMreceiver( ...
                ofdmPHY, ...       %%% 发射机传输的参数
                ofdmPHY.Fs, ...    %   采样
                6*ofdmPHY.Fs, ...  % 上采样
                ofdmPHY.nPkts, ...            % 信道训练长度
                1:1:ofdmPHY.nModCarriers, ...    %导频位置
                g, ...             % 选取一段信号
                ref_seq, ...       % 参考序列
                qam_signal, ...    % qam 矩阵
                'off', ...         % 是否采用CPE
                'off', ...         % 对所有载波进行相位补偿
                'KK');             % 接收方式
            Y2=Signal1(Receiver.ofdmPHY.len*(Receiver.Nr.squ_num-1)+1:Receiver.ofdmPHY.len*Receiver.Nr.squ_num);
            [signal_ofdm_martix_RE1,data_ofdm_martix_RE1,Hf_RE1,data_qam_RE1,qam_bit_RE1]=Receiver.Demodulation(Y2);
            [ber,num]=Receiver.Direcct_Cal_BER(data_qam_RE1(:));

            % 选取性能较好段，进行重新调制
            [ofdm_signal2,~] = nn.ofdm(data_ofdm_martix_RE1);
            % 补上直流
            Re_Signal2=Dc+ofdm_signal2;

            rr2=Y2.'-Re_Signal2;
            % 计算误差
            Error2=cat(1,Error2,rr2);

        end
        Signal2=Signal1-Error2.';
        % 更正解码参数
        Receiver.Total_Preprocessed_signal(ipd_btb);
        %解码
        [ber_total3(index),num_total3(index)]=Receiver.Cal_BER(Signal2);



        % 信号复制
        signal_Re=repmat(Re_Signal,k,1);

        % 残留噪声
        Re=ReceivedSignal-signal_Re.';
        % 更正解码参数
        Receiver.Total_Preprocessed_signal(ipd_btb);
        % 解码(接收信号减去残留噪声)
        [ber_total4(index),num_total4(index)]=Receiver.Cal_BER(ReceivedSignal-Re);


        WB.updata(index);

    end
    WB.closeWaitBar();% 分段解码


end
% % Y=ReceivedSignal-Re;
% % X=Y(Receiver.ofdmPHY.len*(Receiver.Nr.squ_num-1)+1:Receiver.ofdmPHY.len*Receiver.Nr.squ_num);



berplot = BERPlot_David();
% 间隔
berplot.interval=2;
% 字号
berplot.Config.FontSize = 14;
berplot.flagThreshold=1;
berplot.flagRedraw=0;
berplot.flagAddLegend=1;
BER=[ber_total.';ber_total1;ber_total2;ber_total3;ber_total4];
LengendArrary=["80km w/o ","80km w SIC","80km w SIC-2","80km w SIC-3","80km Non-iterative"];
berplot.multiplot(Eb_N0_dB,BER,LengendArrary);
