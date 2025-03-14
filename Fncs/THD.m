% 计算THD

N =  length(signal);
if mod(N,2)
    N = N-1;
    signal(end) = [];
end
freq = fs/N.*[0:N/2-1,-N/2:-1].';
Efout = fftshift(fft(signal));

% If=Efout .* conj(Efout);
% If=pnorm(If);
% IfdBm = 10*log10(If);

% 直流的位置
[~,position_DC]=max(Efout);
% 排除直流后的最大值
[~,position]=max(Efout(position_DC+1:end));
% 信号的频点位置
sig_position=position+position_DC;
% 信号功率
signal_power=abs(Efout(sig_position).^2);
% 将信号幅度设置为0
Y=Efout;
Y(sig_position)=0;
% 杂扰功率
Harmonic_power=sum(abs(Y(position_DC+1:1:end)).^2);
% Cal_THD dB单位
THD_dB = 10*log10(Harmonic_power./signal_power);

% 多音信号的峰值位置查找
peak_att=10;
[pks,locs]=findpeaks(Efout(position_DC+1:end),freq(position_DC+1:end),'MinPeakHeight',peak_att);
Y1=Efout;
Y1(locs+position_DC)=0;
% 计算THD