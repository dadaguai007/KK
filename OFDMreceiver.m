classdef OFDMreceiver < handle

    % 参数直接从 OFDM信号 生成类中 直接引入
    properties
        ofdmPHY; % 传入的信号 参数
        Nr%无量纲参数
        Implementation;% 参考信号实施参数
        Button; % 开关讯号
    end


    methods
        function obj = OFDMreceiver(varargin)
            if numel(varargin) == 11
                obj.ofdmPHY                     = varargin{1} ;% 传递而来的OFDM 参数
                obj.Nr.fOsc                     = varargin{2};
                obj.Nr.fUp                      = varargin{3};
                obj.Nr.nTrainSym                = varargin{4}; % 训练序列长度
                obj.Nr.pilotIndex               = varargin{5};% 导频位置
                obj.Nr.squ_num                  = varargin{6};% 选取第 x 段信号
                obj.Implementation.ref          = varargin{7}; % 参考序列
                obj.Implementation.qam_signal   = varargin{8};% 调制信号参考矩阵

                obj.Button.CPE_Status           = varargin{9};% 默认 关闭 CPE
                obj.Button.PN_Total_Carrier     = varargin{10};% 默认 关闭 所有载波相除相噪
                obj.Button.receive_type         = varargin{11};% 默认 直接接收
                
            else
                error('Number of input variables must be either 0 (default values) or 10');
            end
            obj.Nr.k=1;% 默认为1 ，不变， 选取全部信号时，在函数中更改
        end


        function [compSig,lnSig]=KK_receiver(obj,rxSig)
            % 采样倍数
            m=obj.Nr.fUp/obj.Nr.fOsc;
            N=length(rxSig);
            X=fft(rxSig);
            X_up=[X(1:N/2).'  zeros(1,N*m-N)  X(N/2+1:end).'];
            x_up=m*ifft((X_up),m*N);
            upRxSig=x_up;
            %    upRxSig1=resample(rxSig,fUp,fOsc);
            hn = sqrt(upRxSig);
            % ln of abs
            lnSig = log(abs(hn));

            % phi
            H_sig=hilbert(lnSig);

            % exp
            phi = exp(1i.*(imag(H_sig)));
            % phase recovered signal
            compSig = hn.*phi;
        end
        % 信号预处理
        function  [ReceivedSignal,Dc]=Preprocessed_signal(obj,rxsig)
            if strcmp(obj.Button.receive_type,'KK')
                c=0;
                Dc=mean(rxsig);
                % KK
                [SSB_Sig_KK,~]=obj.KK_receiver(rxsig.'+c);

                %下采样
                ReceivedSignal = downsample(SSB_Sig_KK,obj.Nr.fUp/obj.Nr.fOsc);

            elseif strcmp(obj.Button.receive_type,'DD')
                % 不使用KK算法，使用带宽隔开,需要去除DC
                % DC-remove
                Dc=mean(rxsig);
                rxsig=rxsig-mean(rxsig);
                ReceivedSignal=pnorm(rxsig);
            end

            % 选取某段的信号
            ReceivedSignal=ReceivedSignal(obj.ofdmPHY.len*(obj.Nr.squ_num-1)+1:obj.ofdmPHY.len*obj.Nr.squ_num);
            
        end
        % 信号预处理 （信号全部进行选择）
        function  [ReceivedSignal,Dc]=Total_Preprocessed_signal(obj,rxsig)
            if strcmp(obj.Button.receive_type,'KK')
                c=0;
                Dc=mean(rxsig);
                % KK
                [SSB_Sig_KK,~]=obj.KK_receiver(rxsig.'+c);

                %下采样
                ReceivedSignal = downsample(SSB_Sig_KK,obj.Nr.fUp/obj.Nr.fOsc);

            elseif strcmp(obj.Button.receive_type,'DD')
                % 不使用KK算法，使用带宽隔开,需要去除DC
                % DC-remove
                Dc=mean(rxsig);
                rxsig=rxsig-mean(rxsig);
                ReceivedSignal=pnorm(rxsig);
            end

            % 选取全部信号
            % ReceivedSignal=ReceivedSignal(obj.ofdmPHY.len*(obj.Nr.squ_num-1)+1:obj.ofdmPHY.len*obj.Nr.squ_num);
            obj.Nr.k=length(ReceivedSignal)./obj.ofdmPHY.len;
        end

        function [ber,num]=Cal_BER(obj,ReceivedSignal)
            % 解调算法
            [~,~,~,~,qam_bit]=obj.Demodulation(ReceivedSignal);
            % label信号
            ref_seq =qamdemod(obj.Implementation.ref ,obj.ofdmPHY.M,'OutputType','bit','UnitAveragePower',1);
            ref_seq=ref_seq(:);
            % 计算误码率
            [ber,num,~] = CalcBER(qam_bit,ref_seq);
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        % 输入QAM信号进行解码
        function [ber,num]=Direcct_Cal_BER(obj,ReceivedSignal)
            % 信号解码
            qam_bit=obj.hard_decision(ReceivedSignal);
            % label信号
            ref_seq =qamdemod(obj.Implementation.ref ,obj.ofdmPHY.M,'OutputType','bit','UnitAveragePower',1);
            ref_seq=ref_seq(:);
            % 计算误码率
            [ber,num,~] = CalcBER(qam_bit,ref_seq);
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        % 输入QAM信号进行解码
        function [ber,num]=BER(obj,ReceivedSignal,ref)
            % 信号解码
            qam_bit=obj.hard_decision(ReceivedSignal);
            % label信号
            ref_seq =qamdemod(ref ,obj.ofdmPHY.M,'OutputType','bit','UnitAveragePower',1);
            ref_seq=ref_seq(:);
            % 计算误码率
            [ber,num,~] = CalcBER(qam_bit,ref_seq);
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end
        % 解调
        function [signal_ofdm_martix,data_ofdm_martix,Hf,data_qam,qam_bit]=Demodulation(obj,ReceivedSignal)

            % 解OFDM
            signal_ofdm = reshape(ReceivedSignal, obj.ofdmPHY.fft_size+ obj.ofdmPHY.nCP,[]); % 转换为矩阵形式
            signal_ofdm(1: obj.ofdmPHY.nCP,:) = [];% 去除CP
            signal_ofdm = fft(signal_ofdm);
            % 存储矩阵形式的OFDM 信号
            signal_ofdm_martix=signal_ofdm;
            % get the modulated data carriers
            data_ofdm = signal_ofdm(obj.ofdmPHY.dataCarrierIndex,:);
            % 信道均衡
            [data_ofdm,Hf]=obj.one_tap_equalization(data_ofdm);

            % CPE compensation
            if strcmp(obj.Button.CPE_Status,'on')
                data_ofdm=obj.CPE_Eliminate(data_ofdm);
            end

            if strcmp(obj.Button.PN_Total_Carrier,'on')
                data_ofdm=obj.PN_Total_Eliminate(data_ofdm);
            end

            %保留信号矩阵
            data_ofdm_martix=data_ofdm;
            %归一化
            data_ofdm=data_ofdm(:);% 矩阵转换为行向量
%             data_ofdm = data_ofdm./sqrt(mean(abs(data_ofdm(:)).^2));
            % 硬判决 为 最近的星座点
            data_qam=hard_decision_qam(obj.ofdmPHY.M,data_ofdm);
            % 硬判决 为bit
            qam_bit=obj.hard_decision(data_ofdm);

        end
        
        % 重新生成满足条件的复数信号，  不是平方接收的实数信号
        function ofdm_signal=Remodulation(obj,ReceivedSignal,Dc)
            % 解调获得 信号  和 信道响应
            [~,~,Hf,data_qam,~]=obj.Demodulation(ReceivedSignal);
            % 转换为矩阵形式
            data_qam_martix=reshape(data_qam,obj.ofdmPHY.nModCarriers,[]);
            % 信道响应 叠加
            H_data_qam_martix=data_qam_martix.*repmat(Hf,1,obj.ofdmPHY.nPkts);
            % 重新调制
            % nOffsetSub 行置零 ,positive 放置qam ， 后续置零
            X= ([zeros(obj.ofdmPHY.nOffsetSub,obj.ofdmPHY.nPkts);...
                H_data_qam_martix; ...
                zeros(obj.ofdmPHY.fft_size-obj.ofdmPHY.nModCarriers-obj.ofdmPHY.nOffsetSub,obj.ofdmPHY.nPkts)]);
            % 转换为时域
            ofdmSig=ifft(X);
            % 添加CP
            ofdmSig = [ofdmSig(end-obj.ofdmPHY.nCP+1:end,:);ofdmSig];
            % 并串转换
            ofdmSig = ofdmSig(:);
            % 归一化
%             scale_factor = max(max(abs(real(ofdmSig))),max(abs(imag(ofdmSig))));
%             ofdm_signal = ofdmSig./scale_factor;
            % 还需添加直流项
            ofdm_signal=ofdmSig+Dc;
        end
        
        % 拆分分组
        function [DataGroup,processedGroups,CyclicsGroups]=GroupDemodulation(obj,ReceivedSignal,Grop_index)
            
            % 分组 解码
            DataGroup = cell(1, obj.ofdmPHY.L);
            processedGroups=cell(1,obj.ofdmPHY.L);
            CyclicsGroups=cell(1,obj.ofdmPHY.L);
            for i=1:obj.ofdmPHY.L
                carrier_index=Grop_index(i,:);
                data_group=ReceivedSignal(carrier_index,:);
                DataGroup{i} = data_group;
                if strcmp(obj.Button.Cyclic,'CP_CS')
                    % 去除前缀和后缀载波
                    processedGroups{i}= data_group(obj.ofdmPHY.L_cp+1:end-obj.ofdmPHY.L_cs,:); % 去除CP/CS
                    % 前缀和后缀进行存储
                    CyclicsGroups{i}=cat(1,data_group(1:obj.ofdmPHY.L_cp,:),data_group(end-obj.ofdmPHY.L_cs+1:end,:));
                end
            end
        end

        % 去除CP
        function CP_remove_sig=Remove_CP(obj,ReceivedSignal)
            % 解OFDM
            signal_ofdm = reshape(ReceivedSignal, obj.ofdmPHY.fft_size+ obj.ofdmPHY.nCP,[]); % 转换为矩阵形式
            signal_ofdm(1: obj.ofdmPHY.nCP,:) = [];% 去除CP
            %  矩阵转换为行向量
            CP_remove_sig=signal_ofdm(:);
        end

        % ZF信道均衡
        function [data_kk,Hf]=one_tap_equalization(obj,data_ofdm)
            % channel estimation
            rxTrainSymbol = data_ofdm(:,1:1:obj.Nr.nTrainSym);
            qam_signal_mat=repmat(obj.Implementation.qam_signal,1,1);
            refTrainSymbol = qam_signal_mat(:,1:1:obj.Nr.nTrainSym);
%             % 信道响应
            Hf = mean(rxTrainSymbol./refTrainSymbol,2);

%             % channel equalization
            data_kk = data_ofdm.*repmat(1./Hf,1,obj.ofdmPHY.nPkts*obj.Nr.k);
            %
%             Hf = rxTrainSymbol./refTrainSymbol;
%             data_kk = data_ofdm.*(1./Hf);
        end
        % 相位均衡 CPE
        function data=CPE_Eliminate(obj,data_ofdm)

            phi_mean=angle(mean(data_ofdm(obj.Nr.pilotIndex,:)./...
                obj.Implementation.qam_signal(obj.Nr.pilotIndex,:),1));


            data=data_ofdm.*...
                repmat(exp(-1j.*phi_mean),size(data_ofdm,1),1);
        end

        function data=PN_Total_Eliminate(obj,data_ofdm)

            phi_mean=angle((data_ofdm(obj.Nr.pilotIndex,:)./...
                 obj.Implementation.qam_signal(obj.Nr.pilotIndex,:)));


            data=data_ofdm.*...
                (exp(-1j.*phi_mean));
        end

        % 时域去除PN
        function  [phi_est,data]=Time_Phase_Eliminate(obj,ReceivedSignal,ReconstructSignal)
            % 该方法，不适用直接接收的OFDM系统，接收信号不是所谓的复数信号，而是实数信号，出现较大的误差
            % 去除 CP
            CP_remove_Receiver=obj.Remove_CP(ReceivedSignal);
            CP_remove_Reconstruct=obj.Remove_CP(ReconstructSignal);
            % 重构信号 / 接收信号
            % 相位差估计
            phi_est=angle(CP_remove_Receiver./CP_remove_Reconstruct);             
%             phi_est = angle(conj(CP_remove_Reconstruct) .* CP_remove_Receiver);
%             phi_est=LPF(phi_est,obj.ofdmPHY.Fs,10e9);
            % 补偿
            data=CP_remove_Receiver.*exp(-1j.*phi_est);
            % 转换为矩阵形式
            data_martix= reshape(data, obj.ofdmPHY.fft_size,[]); % 转换为矩阵形式
            % 添加CP
            data_martix = [data_martix(end-obj.ofdmPHY.nCP+1:end,:);data_martix];
            % 并串转换
            data=data_martix(:);
            % 关闭所有消除
            obj.Button.PN_Total_Carrier='off';

        end
         % DD 平方接收 
        function ofdm_time_signal=DD_Type_TimeSignal(obj,signal_ofdm_martix,Dc)
            % 使用 去除CP 后 的频域 OFDM 矩阵
            data_ofdm = signal_ofdm_martix(obj.ofdmPHY.dataCarrierIndex,:);
            X= ([zeros(obj.ofdmPHY.nOffsetSub,obj.ofdmPHY.nPkts);...
                data_ofdm; ...
                zeros(obj.ofdmPHY.fft_size-obj.ofdmPHY.nModCarriers-obj.ofdmPHY.nOffsetSub,obj.ofdmPHY.nPkts)]);
            % 转换为时域
            ofdmSig=ifft(X);
            % 添加CP
            ofdmSig = [ofdmSig(end-obj.ofdmPHY.nCP+1:end,:);ofdmSig];
            % 并串转换
            ofdmSig = ofdmSig(:);
            % 归一化
            scale_factor = max(max(abs(real(ofdmSig))),max(abs(imag(ofdmSig))));
            ofdm_signal = ofdmSig./scale_factor;
            % 还需添加直流项
            ofdm_time_signal=ofdm_signal+Dc;
        end

        function [ber,num]=PCP_Cal_BER(obj,ReceivedSignal,index_data,index_pcp)
            % 解调算法
            [~,~,~,~,qam_bit]=obj.PCP_Demodulation(ReceivedSignal,index_data,index_pcp);
            % label信号
            ref_seq =qamdemod(obj.Implementation.ref ,obj.ofdmPHY.M,'OutputType','bit','UnitAveragePower',1);
            ref_seq=ref_seq(:);
            % 计算误码率
            [ber,num,~] = CalcBER(qam_bit,ref_seq);
            fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);
        end

        % PCP解调
        function [signal_ofdm_martix,data_ofdm_martix,Hf,data_qam,qam_bit]=PCP_Demodulation(obj,ReceivedSignal,index_data,index_pcp)

            % 解OFDM
            signal_ofdm = reshape(ReceivedSignal, obj.ofdmPHY.fft_size+ obj.ofdmPHY.nCP,[]); % 转换为矩阵形式
            signal_ofdm(1: obj.ofdmPHY.nCP,:) = [];% 去除CP
            signal_ofdm = fft(signal_ofdm);
            % 存储矩阵形式的OFDM 信号
            signal_ofdm_martix=signal_ofdm;
            % get the modulated data carriers
            data_ofdm = signal_ofdm(obj.ofdmPHY.dataCarrierIndex,:);
            % 信道均衡
            [data_ofdm,Hf]=obj.one_tap_equalization(data_ofdm);

            % CPE compensation
            if strcmp(obj.Button.CPE_Status,'on')
                data_ofdm=obj.CPE_Eliminate(data_ofdm);
            end

            if strcmp(obj.Button.PN_Total_Carrier,'on')
                data_ofdm=obj.PN_Total_Eliminate(data_ofdm);
            end
            % 相位共轭进行处理
            % 找到索引
            data_index=index_data-obj.ofdmPHY.nOffsetSub;
            pcp_index=index_pcp-obj.ofdmPHY.nOffsetSub;
            % 选取共轭数据
            data_ofdm_data=data_ofdm(data_index,:);
            data_ofdm_pcp=data_ofdm(pcp_index,:);

            %   共轭消除 得到恢复信号
            data_ofdm=(data_ofdm_data+conj(data_ofdm_pcp))/2;
            %保留信号矩阵
            data_ofdm_martix=data_ofdm;
            %归一化
            data_ofdm=data_ofdm(:);% 矩阵转换为行向量
%             data_ofdm = data_ofdm./sqrt(mean(abs(data_ofdm(:)).^2));
            % 硬判决 为 最近的星座点
            data_qam=hard_decision_qam(obj.ofdmPHY.M,data_ofdm);
            % 硬判决 为bit
            qam_bit=obj.hard_decision(data_ofdm);

        end


        % 硬判决
        function data_qam=hard_decision(obj,Receiver_data)
            data_qam = qamdemod(Receiver_data,obj.ofdmPHY.M,'OutputType','bit','UnitAveragePower',1);
            data_qam=data_qam(:);
        end
        % EVM测算

        function [rmsEVM_symbol,rmsEVM_subcarrier,rmsEVM_martix]=EVM_Measure_martix(obj,data_ofdm_martix)
            %EVM measure
            evm = comm.EVM(AveragingDimensions=1);
            rmsEVM_symbol = evm(data_ofdm_martix,obj.Implementation.qam_signal);


            evm = comm.EVM(AveragingDimensions=2);
            rmsEVM_subcarrier = evm(data_ofdm_martix,obj.Implementation.qam_signal);


            evm = comm.EVM(AveragingDimensions=[1 2]);
            rmsEVM_martix = evm(data_ofdm_martix,obj.Implementation.qam_signal);

            fprintf('EVM = %1.7f\n',rmsEVM_martix);
        end

        function rmsEVM_subcarrier=EVM_Measure_symbol(obj,data_ofdm_martix,symbol_indedx)
            %EVM measure
            % 某个符号上所有载波的EVM
            evm = comm.EVM(AveragingDimensions=2);
            rmsEVM_subcarrier = evm(data_ofdm_martix(:,symbol_indedx),obj.Implementation.qam_signal(:,symbol_indedx));

        end

        % 星座图绘制
        function scatter_plot(obj,ReceivedSignal)
            [~,data_ofdm_martix,~,~,~]=obj.Demodulation(ReceivedSignal);
            dara_ofdm_squ=data_ofdm_martix(:);
            scatterplot(dara_ofdm_squ);
        end

    end


end