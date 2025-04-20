addpath('Output\');
addpath('Plot\')
load('BER_Dither_0_2_5_8_10.mat')
load('BER_Algriom_Dither_decode_Group_0_2_5_8_10.mat')
load('BER_Algriom_Dither_decode_one_by_one_0_2_5_8_10.mat')
EbNo=15:30;
berplot = BERPlot_David();
% 间隔
berplot.interval=3;
% 字号
berplot.Config.FontSize = 14;
berplot.flagThreshold=1;
berplot.flagRedraw=0;
berplot.flagAddLegend=1;

BER=[ber_total(1,:);...
    ber_total(2,:);ber_total_iter(2,:);ber_total_iter_Group(2,:);];
LengendArrary=["80km w/o 0%Vpi",...
    "80km w/o 2%Vpi","80km w SIC Total 2%Vpi","80km w SIC Group 2%Vpi",];
berplot.multiplot(EbNo,BER,LengendArrary);

% berplot.plot(EbNo,ber_total(1,:));