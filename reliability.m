% 置信度表示

function gamma=reliability(y)

% 星座图数量
M=4;
% 星座图 ，归一化
const = real(pnorm(pammod(0:M-1,M)));

gamma=zeros(length(y),1);
for index=1:length(y)
    input_y = y(index);

    if abs(input_y)<max(const)
        % 进行硬判决
        y_hat=hard_decision_pam(M,input_y);
        gamma(index)=1-abs(input_y-y_hat);
    else
        gamma(index)=1;
    end

end

end
