function [compSig,lnSig] = KK(rxSig,fOsc,fUp)
   % up-sampling
   upRxSig = resample(rxSig,fUp,fOsc); % 5 sps

% m=fUp/fOsc;
% N=length(rxSig);
% X=fft(rxSig);
% X_up=[X(1:N/2).'  zeros(1,N*m-N)  X(N/2+1:end).'];
% x_up=m*ifft((X_up),m*N);
% upRxSig=x_up.';

% square root
   hn = sqrt(upRxSig);
   % ln of abs
   lnSig = log(abs(hn));
   % fft
   fftSig = fft(lnSig);
   % i*sign(w)
   N = length(fftSig);
   freq = [0:N/2-1,-N/2:-1].';
   signFFTSig = -1i.*sign(freq).*fftSig;
   % ifft
   ifftSig = ifft(signFFTSig);

   % exp (此时应该选取信号经过希尔伯特变化后的实部)
   phi = exp(1i.*(real(ifftSig)));
   % phase recovered signal
   compSig = hn.*phi;
end