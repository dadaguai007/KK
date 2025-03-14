function h=Window_FIR_Low_High_Desgin(fs,fc,fir_type,window_type,flag_plot)

if nargin<4
    flag_plot=1;
end

% fs=8000;         %采样频率
% fc=[1000 1500];  %过渡带   低通：前一个为截止项；高通：后一个为起始项

if strcmp(fir_type,'low')
    mag=[1 0];       %窗函数的理想滤波器幅度
    dev=[0.01 0.05]; %纹波
    [n,wn,beta,ftype]=kaiserord(fc,mag,dev,fs);  %获取凯塞窗参数

elseif strcmp(fir_type,'high')
    mag=[0 1];       %窗函数的理想滤波器幅度
    dev=[0.01 0.05]; %纹波
    [n,wn,beta,ftype]=kaiserord(fc,mag,dev,fs);  %获取凯塞窗参数
end

if strcmp(window_type,'kaiser')

    h_kaiser=fir1(n,wn,ftype,kaiser(n+1,beta));
    %求滤波器的幅频响应
    m_kaiser=20*log(abs(fft(h_kaiser,1024)))/log(10);
    filterResponse=m_kaiser;
    h=h_kaiser;
elseif strcmp(window_type,'hanmm')
    if strcmp(fir_type,'low')
        h_hamm=fir1(n,fc(2)*2/fs);
    else
        h_hamm=fir1(n,fc(2)*2/fs,'high');
    end
    %求滤波器的幅频响应
    m_hamm=20*log(abs(fft(h_hamm,1024)))/log(10);
    filterResponse=m_hamm;
    h=h_hamm;
elseif strcmp(window_type,'optimal')
    if strcmp(fir_type,'low')
        fpm=[0 fc(1)*2/fs fc(2)*2/fs 1];  %firpm函数的频段向量
        magpm=[1 1 0 0];                  %firpm函数的幅值向量
    else
        fpm=[0 fc(1)*2/fs fc(2)*2/fs 1];  %firpm函数的频段向量
        magpm=[0 0 1 1];                  %firpm函数的幅值向量
    end
    %设计最优滤波器
    h_pm=firpm(n,fpm,magpm);
    %求滤波器的幅频响应
    m_pm=20*log(abs(fft(h_pm,1024)))/log(10);
    filterResponse=m_pm;
    h=h_pm;
end

%设置幅频响应的横坐标单位为Hz
x_f= fs * (-0.5:1/length(m_kaiser):0.5-1/length(m_kaiser));


if flag_plot
    figure;
    plot(x_f./1e9,filterResponse,'b');
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    box on;
    set(gca, 'FontName', 'Arial', 'FontSize', 14);
    set(gca, 'LineWidth', 1.25);
end
% filter(channel，1，x)滤波器的使用：没有延时效果 
end

