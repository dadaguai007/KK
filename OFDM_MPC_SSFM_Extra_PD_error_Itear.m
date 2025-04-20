% dither 对 多载波信号的影响

% 迭代消除

clear;close all;clc;
addpath('Fncs\')
% addpath('D:\PhD\Project\Base_Code\Base\')
addpath('D:\BIT_PhD\Base_Code\Codebase_using\')

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
    ipd_btb=pnorm(ipd_btb);
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

    % 信号预处理
    [ReceivedSignal,dc]=Receiver.Total_Preprocessed_signal(ipd_btb);
    % 场信号的载波
    Dc=mean(ReceivedSignal);

    % BER 计算
    [ber_total,num_total]=Receiver.Cal_BER(ReceivedSignal);


    % 设置单边带信号
    VbQ1 = Vdither(1)*Creat_dither1(Fs_new,f2,N*(f2/f1));
    VbQ_ssb = 1j*Vdither(1)*Creat_ssb(Fs_new,f2,N*(f2/f1));
    VbI_ssb = Vdither(1)*Creat_ssb(Fs_new,f1,N);
    %     VbQ1 = 0.1*Creat_dither1(Fs_new,f2,N*(f2/f1));
    %     VbQ_ssb = 1j*0.1*Creat_ssb(Fs_new,f2,N*(f2/f1));
    %     VbI_ssb = 0.1*Creat_ssb(Fs_new,f1,N);
    ipd_pd=ipd_btb;
    for jj=1:60
        dd=real(Dc);
        % 重新计算E2和E1
        Rece_remove_dc=ReceivedSignal-dd;

        % 载波与dither 拍频
        E1=2*real(dd)*VbI+2*real(dd)*VbQ1;
        % 信号与dither拍频
        I_beat=Rece_remove_dc.*conj(VbI_ssb)+conj(Rece_remove_dc).*VbI_ssb;
        Q_beat=Rece_remove_dc.*conj(VbQ_ssb)+conj(Rece_remove_dc).*VbQ_ssb;

        %     E1=2*real(dd)*(VbI);
        E2=I_beat+Q_beat;

        alpha=0.002;
        ipd_error=alpha*(E1+E2);

        ipd_pd=ipd_pd-ipd_error;


        % 信号预处理
        [ReceivedSignal,~]=Receiver.Total_Preprocessed_signal(ipd_pd);
        Dc=mean(ReceivedSignal);

        % BER 计算
        [ber_total1,num_total1]=Receiver.Cal_BER(ReceivedSignal);
    end

% 选取局部信号进行相应的操作

    % 发射机参数
    ofdmPHY=nn;
    for i=5

        %%---------------------------------------        解码       ---------------------------%%
        Receiver=OFDMreceiver( ...
            ofdmPHY, ...       %%% 发射机传输的参数
            ofdmPHY.Fs, ...    %   采样
            6*ofdmPHY.Fs, ...  % 上采样
            10, ...            % 信道训练长度
            1:1:ofdmPHY.nModCarriers, ...    %导频位置
            i, ...             % 选取第一段信号
            ref_seq, ...       % 参考序列
            qam_signal, ...    % qam 矩阵
            'off', ...         % 是否采用CPE
            'off', ...         % 对所有载波进行相位补偿
            'KK');             % 接收方式

        % 信号预处理
        [receive,~]=Receiver.Preprocessed_signal(ipd_btb);

        % BER 计算
        [ber,num]=Receiver.Cal_BER(receive);
        % 选取某段的信号
        rx_pd=ipd_btb(Receiver.ofdmPHY.len*(Receiver.Nr.squ_num-1)+1:Receiver.ofdmPHY.len*Receiver.Nr.squ_num);
        % 选取某一段
        rx_VbI=VbI(Receiver.ofdmPHY.len*(Receiver.Nr.squ_num-1)+1:Receiver.ofdmPHY.len*Receiver.Nr.squ_num);
        rx_VbQ1=VbQ1(Receiver.ofdmPHY.len*(Receiver.Nr.squ_num-1)+1:Receiver.ofdmPHY.len*Receiver.Nr.squ_num);
        rx_VbI_ssb=VbI_ssb(Receiver.ofdmPHY.len*(Receiver.Nr.squ_num-1)+1:Receiver.ofdmPHY.len*Receiver.Nr.squ_num);
        rx_VbQ_ssb=VbQ_ssb(Receiver.ofdmPHY.len*(Receiver.Nr.squ_num-1)+1:Receiver.ofdmPHY.len*Receiver.Nr.squ_num);
        for kk=1:20
            dd=real(Dc);
            % 重新计算E2和E1
            rece_remove_dc=receive-dd;

            % 载波与dither 拍频
            rx_E1=2*real(dd)*rx_VbI+2*real(dd)*rx_VbQ1;
            % 信号与dither拍频
            rx_I_beat=rece_remove_dc.*conj(rx_VbI_ssb)+conj(rece_remove_dc).*rx_VbI_ssb;
            rx_Q_beat=rece_remove_dc.*conj(rx_VbQ_ssb)+conj(rece_remove_dc).*rx_VbQ_ssb;

            %     E1=2*real(dd)*(VbI);
            rx_E2=rx_I_beat+rx_Q_beat;

            alpha=0.002;
            rx_error=alpha*(rx_E1+rx_E2);

            rx_pd=rx_pd-rx_error;


            % 信号预处理
            [receive,dc2]=Receiver.Total_Preprocessed_signal(rx_pd);
            Dc=mean(receive);
            % BER 计算
            [ber_total2,num_total2]=Receiver.Cal_BER(receive);
        end


    end

end



[If_pd,fre]=mon_ESA_flag(ipd_btb,fs,0);
[If_dbm_error,fre]=mon_ESA_flag(ipd_error,fs,0);
[If_dbm,fre]=mon_ESA_flag(ipd_pd,fs,0);
[If_kk,fre]=mon_ESA_flag(ReceivedSignal,fs,0);

figure;hold on;

plot(fre,If_kk)
plot(fre,If_pd)
plot(fre,If_dbm_error)
plot(fre,If_dbm)
legend('恢复算法','接收','误差','接收减去误差')
xlim([11.5815 11.5835]);