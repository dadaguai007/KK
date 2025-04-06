% Volterra DFE LMS
function [y,e,w]=ReLU_dfe_lms(xn,dn,sps,ref,taps_list,step_len)
%三阶volterra
% 转换为列
xn = xn(:);
dn = dn(:);
% 序列长度
L1 = length(xn);
% 期望输入
L2 = length(dn);
% 信号输出长度
n = 2*floor(L1/sps/2);
% 输出
y = zeros(n, 1);


% 前馈输入的抽头
% 一阶
tapslen_1=taps_list(1);
% 二阶
tapslen_2=taps_list(2);
% 三阶
tapslen_3=taps_list(3);
%反馈输入抽头
fblen_1=taps_list(4);
fblen_2=taps_list(5);
fblen_3=taps_list(6);

% 每个长度都是不同的：整个数组长度，从数组抽取两个元素的长度，从数组中抽取三个元素的长度

% 抽头长度
w = zeros(tapslen_1+tapslen_2 + tapslen_3 +...
    fblen_1+fblen_2 + fblen_3 ,1);


% 参考抽头 % 参考抽头应该存在两种方法进行定义
xn = cat(1,xn(ref+1:end),dn(end-ref+1:end));
% 反馈抽头
fb = zeros(1,fblen_1);
% 前馈抽头
x1=zeros(tapslen_1,1);

for idx = 1:n-1
    %构建volterra输入 
    %一阶前馈输入
    x1 = cat(1,x1(sps+1:end),xn(sps*idx-sps+1:1:sps*idx));
    %二阶前馈输入 % 从x1中确定好抽头的数量
    x2 = x1(round((tapslen_1-tapslen_2)/2)+1 : end - fix((tapslen_1-tapslen_2)/2));
%     x2 = x1(1:tapslen_2);
    % 二阶核
    x2_vol=Build_ReLU_Input(x2);

    %三阶前馈输入
    x3 = x1(round((tapslen_1-tapslen_3)/2)+1 : end - fix((tapslen_1-tapslen_3)/2));
%     x3=x1(1:tapslen_3);
    % 三阶核
    x3_vol=Build_ReLU_Input(x3);



    % 反馈输入
    %一阶反馈输入
    fb1_vol = fb(1:fblen_1);
    %二阶反馈输入 

    fb2_vol=Build_ReLU_Input(fb(1:fblen_2));
    %三阶反馈输入
    fb3_vol=Build_ReLU_Input(fb(1:fblen_3));
    %组合所有输入
    x=x1.';
    
    x_all = [x x2_vol x3_vol fb1_vol fb2_vol fb3_vol];
    y(idx)=x_all * w;
    if idx+tapslen_1-1>L2
        % DFE 的判决
        Y_D=DFE_decision(y(idx));
        y_d(idx)=Y_D;
        fb_new = Y_D;
        e(idx) = Y_D - y(idx);
    else
        e(idx) = dn(idx) - y(idx);
        fb_new = dn(idx);
    end

    % 抽头更新
    %使用lms更新抽头
    w = w +  e(idx) * step_len * x_all.';
    %反馈更新
    fb = [fb(2:end) fb_new];

end

% 进行聚类（好像是直接进行聚类）


end