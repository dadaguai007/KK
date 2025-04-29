function [ipd_error]=iteraError(E,fs,alpha,Dc,Vdither,f1,f2)

% f1=40e3;
% f2=60e3;
N=length(E)/(fs/f1);

% 创建dither信号
VbQ_sin = Vdither*Creat_dither1(fs,f2,N*(f2/f1));
VbI_sin = Vdither*Creat_dither1(fs,f1,N);
VbI_cos = Vdither*Creat_dither(fs,f1,N);
VbQ_cos = Vdither*Creat_dither(fs,f2,N*(f2/f1));


% 除去直流
E_removedc=E;


% 载波与dither 拍频
E5=real(Dc)*VbI_cos+real(Dc)*VbQ_sin;
E4=real(Dc)*VbI_cos-real(Dc)*VbQ_sin;
E1=E5+E4;



%  负频率 拍频
I_beat=real(E_removedc).*VbI_cos-imag(E_removedc).*VbI_sin;
Q_beat=real(E_removedc).*VbQ_sin+imag(E_removedc).*VbQ_cos;

E2=I_beat+Q_beat;


% 正频率
I_beat1=real(E_removedc).*VbI_cos+imag(E_removedc).*VbI_sin;
Q_beat1=-real(E_removedc).*VbQ_sin+imag(E_removedc).*VbQ_cos;

E3=I_beat1+Q_beat1;


% 自拍频 
%PD
paramPD.B =fs;
paramPD.R =1;
paramPD.type = 'ideal';
paramPD.Fs=fs;

Dither=VbI_cos+1j*VbI_sin;
%pd
E6 = pd(Dither, paramPD);


% 消除不需要的频率分量
ipd_error=alpha*(E1+E2+E3+E6);


end