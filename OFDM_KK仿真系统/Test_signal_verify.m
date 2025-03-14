
addpath('Fncs\')
addpath('D:\001-处理中\相干算法\optical_communication')
clc;clear;close all;
rng(10);
M=16;
SpS=4;
Nbits = (log2(M)*1e6*5);
Rs  = 20e9; 
Fs = Rs*SpS;
psfRollOff=0.01;
psfLength=100;
psfShape='sqrt';
%设计根升余弦脉冲成型滤波器
hsqrt = rcosdesign(psfRollOff,psfLength,SpS,psfShape);
% Generate random bits
bitsTx = randi([0, 1], log2(M), Nbits/log2(M));

% Map bits to constellation symbols
symbTx=qammod(bitsTx,M,'InputType','bit','UnitAveragePower',1);
symbTx = pnorm(symbTx);

% symbTx=awgn(symbTx,20,'measured');
phi=pi/5.5;
symbTx_phase=symbTx.*exp(1j*phi);


[~,t_up]=freq_time_set(Nbits/log2(M),Fs);
f_lo=0.3e4;
symbTx_fre=symbTx.*exp(-1j*(2*pi*f_lo*t_up-phi));
symbTx_fre=awgn(symbTx_fre,30,'measured');
scatterplot(symbTx_fre);
% scatplot(real(symbTx_fre(1:1e6)),imag(symbTx_fre(1:1e6)));
figure;
plot(real(symbTx_fre),'.')