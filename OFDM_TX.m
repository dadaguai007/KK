addpath(genpath('Fcns'));

nn=OFDMQAMN();
nn.DataType='rand';%两个选项：prbs，rand
nn.NSym = 30000;
nn.fft_size = 1024;
nn.nPkts = 100;
nn.nCP = 32;
nn.nModCarriers = 300;
nn.nOffsetSub =3; 
nn.order = 10;
nn.M = 16;
nn.prbsOrder = 11;
nn.Rs = 64e9;
nn.Fs = 64e9;
nn.Nsam =nn.Fs/nn.Rs ;
nn.psfRollOff=0.01;
nn.psfLength=256;
nn.psfShape='sqrt';
nn.psfshape='Raised Cosine';
nn.len= (nn.fft_size+nn.nCP)*nn.nPkts; % OFDM 符号长度
nn.dataCarrierIndex=nn.nOffsetSub+(1:nn.nModCarriers);
f_idex=((nn.nModCarriers+nn.nOffsetSub)/nn.fft_size)*nn.Fs; 

% OFDM：1024*8+32*8