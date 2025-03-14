
figure('Visible', showtag);
subplot(3,1,1)
plot(real(SSB_Sig_KK));hold on;
plot(imag(SSB_Sig_KK));
title('IQ输出')
subplot(3,1,2)
plot(ipd_btb)
title('PD接收')
subplot(3,1,3)
plot(1:idx,number_num);
title('误码数')
xlabel('符号数')
if savefigtag
    S=num2str(snr(index));
    datapath =strcat('Output\','Fig_',S,'dB',pre);
    if ~exist(datapath,'dir')
        mkdir(datapath);
    end
    saveas(gcf,sprintf('%s\\Error%s_%s.png',datapath,pre,vpp));
end

figure('Visible', showtag);
subplot(3,1,1)
plot(rmsEVM_subcarrier)
title('载波分布EVM')
subplot(3,1,2)
plot(rmsEVM_symbol)
title('符号分布EVM')
subplot(3,1,3)
plot(S_l)
xlabel('Index of OFDM symbol')
title('符号的相位变化')

if savefigtag

    S=num2str(snr(index));
    saveas(gcf,sprintf('%s\\EVM%s_%s.png',datapath,pre,vpp));
end