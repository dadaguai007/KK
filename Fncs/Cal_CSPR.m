function Eout1=Cal_CSPR(sigTxo,Ai,Vpi,phi)

% Cal CSPR
param.Vpi=Vpi;
param.VbI=-phi*Vpi;
param.VbQ=-phi*Vpi;
param.Vphi=Vpi/2;
Eout1 = iqm(Ai, 0, param);
Eout_s=sigTxo-Eout1*ones(1,length(sigTxo));
power1=signalpower(Eout1);
Ps=signalpower(Eout_s);
CSPR = power1/(Ps);
CSPR=10*log10(CSPR);
fprintf("the CSPR is %1.2fdB\n",CSPR);

end