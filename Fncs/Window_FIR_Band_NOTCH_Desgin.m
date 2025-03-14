function h=Window_FIR_Band_NOTCH_Desgin(fs,fc,fir_type,window_type,flag_plot)

if nargin<4
    flag_plot=1;
end

% fs=8000;         %采样频率
% fc=[1000 1500 3000 3500];  %过渡带   带通：第二个频率为起始 第三个频率为截止
%                            带阻：第二个频率为起始 第三个频率为截止
if strcmp(fir_type,'band')
    mag=[0 1 0];       %窗函数的理想滤波器幅度
    dev=[0.01 0.05 0.01]; %纹波
    [n,wn,beta,ftype]=kaiserord(fc,mag,dev,fs);  %获取凯塞窗参数
    n = n + rem(n,2);
elseif strcmp(fir_type,'notch')
    mag=[1 0 1];       %窗函数的理想滤波器幅度
    dev=[0.05 0.0001 0.05]; %纹波
    [n,wn,beta,ftype]=kaiserord(fc,mag,dev,fs);  %获取凯塞窗参数
    n = n + rem(n,2);
end

if strcmp(window_type,'kaiser')

    h_kaiser=fir1(n,wn,ftype,kaiser(n+1,beta));
    %求滤波器的幅频响应
    m_kaiser=20*log(abs(fft(h_kaiser,1024)))/log(10);
    filterResponse=m_kaiser;
    h=h_kaiser;
elseif strcmp(window_type,'hanmm')
    if strcmp(fir_type,'band')
        h_hamm=fir1(n,[fc(2)*2/fs fc(3)*2/fs],"DC-0");
    else
        h_hamm=fir1(n,[fc(2)*2/fs fc(3)*2/fs],'stop');
    end
    %求滤波器的幅频响应
    m_hamm=20*log(abs(fft(h_hamm,1024)))/log(10);
    filterResponse=m_hamm;
    h=h_hamm;
elseif strcmp(window_type,'optimal')
    if strcmp(fir_type,'band')
        fpm=[0 fc(1)*2/fs fc(2)*2/fs  fc(3)*2/fs fc(4)*2/fs 1];  %firpm函数的频段向量
        magpm=[0 0 1 1 0 0];                  %firpm函数的幅值向量
    else
        fpm=[0 fc(1)*2/fs fc(2)*2/fs  fc(3)*2/fs fc(4)*2/fs 1];  %firpm函数的频段向量
        magpm=[1 1 0 0 1 1];                  %firpm函数的幅值向量
    end
    %设计最优滤波器
    h_pm=firpm(n,fpm,magpm);
    %求滤波器的幅频响应
    m_pm=20*log(abs(fft(h_pm,1024)))/log(10);
    filterResponse=m_pm;
    h=h_pm;
end


if flag_plot
    [H,f] = freqz(h,1,1024,fs);
    figure;
    plot(f,abs(H))
    grid
end
% filter(channel，1，x)滤波器的使用：没有延时效果 
end

