% CSPR值计算
% 偏离程度
a=phi;
Vpi=paramIQ.Vpi;

paramIQ.Vpi=9;
paramIQ.VbI=-phi*paramIQ.Vpi;
paramIQ.VbQ=-phi*paramIQ.Vpi;
paramIQ.Vphi=paramIQ.Vpi/2;

% theory
Eout1 = iqm(Ai, 0, paramIQ);
Eout_s=sigTxo-Eout1*ones(1,length(sigTxo));
power1=signalpower(Eout1);
Ps=signalpower(Eout_s);
CSPR = power1/(Ps);
CSPR=10*log10(CSPR);
fprintf('CSPR1 = %1.7f\n',CSPR);
U=mean(ipd_btb);
if 0
% Bo-simply CSPR
U=mean(ipd_btb);
V=var(ipd_btb);
x1=sqrt(U.^2-V);
x2=U-x1;
CSPR1 = x1/x2;
CSPR1=10*log10(CSPR1);
fprintf('CSPR3 = %1.7f\n',CSPR1);
fprintf('the index number mean = %1.7f\n',U);

% Sim
c=(pi.^2)/(4*(Vpi.^2));
d=1/2*(1-cos(a*pi));
d1=1/2*(1+cos(a*pi));
p=var(Amp*real(signal));

Mindx = str2double(sprintf('%.3f', m/2));
p1=((Mindx*Vpi).^2)/9;
CSPR4 = d1/(p1*c*d) ;
CSPR4=10*log10(CSPR4);
fprintf('CSPR4 is %1.3f\n',CSPR4);
end