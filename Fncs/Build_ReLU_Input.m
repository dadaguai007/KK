function [res] = Build_ReLU_Input(input_X)
% 构建 ReLU-based 输入
% input_X: 输入信号向量
% order: Volterra 级数
res = [];
u=1;


for i=1:length(input_X)
    nonlinear=max(0,input_X(i));
    res=[res nonlinear];
    u=u+1;
end

end
