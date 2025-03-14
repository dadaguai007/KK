addpath(genpath('Fcns'));

nn=OFDMQAMN();
nn.DataType='rand';%两个选项：prbs，rand
nn.NSym = 200000;
nn.fft_size = 512;
nn.nPkts = 1000;
nn.nCP = 32;
nn.nModCarriers = 200;
nn.nOffsetSub =3; 
nn.order = 2;
nn.M = 16;
nn.prbsOrder = 11;
nn.Rs = 1e9;
nn.Fs = 1e9;
nn.Nsam =nn.Fs/nn.Rs ;
nn.psfRollOff=0.01;
nn.psfLength=256;
nn.psfShape='sqrt';
nn.psfshape='Raised Cosine';

% OFDM：1024*8+32*8