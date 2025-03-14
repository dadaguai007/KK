function [compSig,lnSig,H_sig] = KK_MPC(rxSig,fOsc,fUp)
% up-sampling (频域上采样)

m=fUp/fOsc;
N=length(rxSig);
X=fft(rxSig);
X_up=[X(1:N/2).'  zeros(1,N*m-N)  X(N/2+1:end).'];
x_up=m*ifft((X_up),m*N);
upRxSig=x_up;
%    upRxSig=rxSig;
hn = sqrt(upRxSig);
% ln of abs
lnSig = log(abs(hn));

% phi
H_sig=hilbert(lnSig);

% exp
phi = exp(1i.*(imag(H_sig)));
% phase recovered signal
compSig = hn.*phi;

end