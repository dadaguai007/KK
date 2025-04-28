function [recoverI,ipd_error]=iteraElimate_ssfm(E,ipd_pd,alpha,Dc,VbQ_sin,VbQ_cos,VbI_sin,VbI_cos)


% 除去直流
E_removedc=E-real(Dc);


% 载波与dither 拍频
E5=real(Dc)*VbI_cos+real(Dc)*VbQ_sin;
E4=real(Dc)*VbI_cos-real(Dc)*VbQ_sin;
E1=E5+E4;



%  负频率 拍频
I_beat=real(E_removedc).*VbI_cos-imag(E_removedc).*VbI_sin;
Q_beat=real(E_removedc).*VbQ_sin+imag(E_removedc).*VbQ_cos;

E2=I_beat+Q_beat;

% 消除不需要的频率分量
ipd_error=alpha*(E1+E2);
% ipd_error=alpha*(E2);
recoverI=ipd_pd-ipd_error;

end