function [out,tout,fout]=quickSono(SUTf,t,f, windLen,overlapPct)



oneWindInd=(1:windLen)-round(windLen/2);
increment=round((1-overlapPct)*windLen);
windMids=round(windLen/2):increment:numel(SUTf)-round(windLen/2)-1; % center of each stft window
inds2d=oneWindInd+windMids';
SUT2D=SUTf(inds2d').*hann(windLen);
dt=mean(diff(t));
out=fftshift(ifft(ifftshift(SUT2D)))/dt;
out=out';
tout=t(windMids);
fout=linspace(-1/(2*dt),1/(2*dt),windLen);


end