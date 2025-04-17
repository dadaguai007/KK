function [x_clipped,x_mean]=clipping_mean(x,CL,x_mean)
% CL   : Clipping Level
% sigma: sqrt(variance of x)
% 可尝试百分之1

if nargin<3
    %求取均值
  x_mean=mean(x); 

  %sigma=std(x)
end
CL = CL*x_mean;
x_clipped = x; 
% 小于规定值，进行削波
ind = find(abs(x)<CL); % Indices to clip
% % x_clipped(ind) = x(ind)./abs(x(ind))*CL;% 考虑了符号

x_clipped(ind) = CL;
