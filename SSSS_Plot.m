clc;clear;close all;
load('D:\PhD\Project\KK算法\Simulate\MPC_KK\Output\BER_Dither_w_o_40k_btb.mat');
BTB=ber_total;
load('D:\PhD\Project\KK算法\Simulate\MPC_KK\Output\BER_Dither_Itera_40k_btb.mat')
ber_sic=ber_total1;
load('D:\PhD\Project\KK算法\Simulate\MPC_KK\Output\BER_Dither_Itera_40k_btb_label.mat')
ber_label=ber_total1;
load('D:\PhD\Project\KK算法\Simulate\MPC_KK\Output\BER_Dither_Itera_40k_btb_re_signal.mat')
ber_sic_remode=ber_total1;
Eb_N0_dB=15:30;

berplot = BERPlot_David();
% 间隔
berplot.interval=3;
% 字号
berplot.Config.FontSize = 14;
berplot.flagThreshold=1;
berplot.flagRedraw=0;
berplot.flagAddLegend=1;
BER=[BTB.';ber_total.';ber_group_total1;ber_sic;ber_label;ber_sic_remode];
LengendArrary=["w/o dither","w dither Total","w dither Group","w dither SIC","w dither SIC using label",...
    "w dither SIC by remod"];
berplot.multiplot(Eb_N0_dB,BER,LengendArrary);