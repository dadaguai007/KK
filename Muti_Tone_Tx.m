% 多音信号生成
addpath('Fncs\')

% 信号幅度
As=1;
% 频率
f_inter=5e6;
SpS=20;
% 采样率
fs=f_inter*SpS;
% Time array
T=1/fs;

N=4e5;
% 创建时间轴
[~,t_up]=freq_time_set(N,fs);
for i=1:6
f_index(i)=f_inter*i;
end
% 角频率
w=2*pi*f_index;
signal=0;
for i=1:6
signal=signal+As*exp(1j*w(i)*t_up);
end

