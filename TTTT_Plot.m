clc;clear;close all;

load('D:\PhD\Project\KK算法\Simulate\MPC_KK\Output\BER_Dither_0_2_5_8_10.mat')
BTB=ber_total(1,:);
ber_total=ber_total(2,:);
load('D:\PhD\Project\KK算法\Simulate\MPC_KK\Output\BER_Dither_Itera_40k_ssfm_re_signal.mat')
ber_sic_remod=ber_total1;
load('D:\PhD\Project\KK算法\Simulate\MPC_KK\Output\BER_Algriom_Dither_0_2_5_8_10.mat')
ber_sic=ber_total_iter(2,:);
load('D:\PhD\Project\KK算法\Simulate\MPC_KK\Output\BER_Algriom_Dither_decode_Group_0_2_5_8_10.mat')
ber_Group=ber_total_iter_Group(2,:);
Eb_N0_dB=15:30;

berplot = BERPlot_David();
% 间隔
berplot.interval=3;
% 字号
berplot.Config.FontSize = 14;
berplot.flagThreshold=1;
berplot.flagRedraw=0;
berplot.flagAddLegend=1;
BER=[BTB;ber_total;ber_Group;ber_sic;ber_sic_remod];
LengendArrary=["80km w/o dither","80km dither Total","80km dither Group","80km dither SIC","80km dither SIC by remod"];
berplot.multiplot(Eb_N0_dB,BER,LengendArrary);