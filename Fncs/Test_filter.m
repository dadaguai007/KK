% 滤波器的使用代码

% N=length(s_recovery)+1;
% fc_stop=[15e3 200e3];
% fir_type='stop';
% FIR_plot=1;
% FIR_use;
% 频域处理 可以用filter进行测试
% dataout1=ifft(fft(s_recovery).*fft(b_stop(1:end-1))) ;

%
% s_recovery=NTF(s_recovery,fs,5e3,800e3,0);


% 数据量大，算不动
% fc_stop=[10e3 200e3];
% fil_type='notch';
% dataout1=Filter_Work(s_recovery,fc_stop,fs,fil_type);

%
% % 算不动，optimal更为算不动
% fir_type='notch';
% window_type='kaiser';
% fc=[10e3 40e3 60e3 70e3];
% h=Window_FIR_Band_NOTCH_Desgin(fs,fc,fir_type,window_type,1);
% dataout1=filter(h,1,s_recovery);