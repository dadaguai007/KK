% 置信度-非线性函数构造 
function f=reliability_nonlinear(x,a,b)

f_up=1-exp(-a*(x/b-1));
f_down=1+exp(-a*(x/b-1));

f=1/2*(f_up./f_down+1);

end