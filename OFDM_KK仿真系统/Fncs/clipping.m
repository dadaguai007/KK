function [x_clipped,sigma]=clipping(x,CL,sigma)
% CL   : Clipping Level
% sigma: sqrt(variance of x)

if nargin<3
  x_mean=mean(x); 
  x_dev=x-x_mean; 
  sigma=sqrt(x_dev*x_dev'/length(x));
end
CL = CL*sigma;
x_clipped = x;  
ind = find(abs(x)>CL); % Indices to clip
x_clipped(ind) = x(ind)./abs(x(ind))*CL;
