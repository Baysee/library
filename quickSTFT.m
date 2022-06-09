function [out,tout,fout]=quickSTFT(SUT,t,f, windLen,overlapPct)


% Convert windlen to indices
windLen=round(windLen/mean(diff(t)));
oneWindInd=(1:windLen)-round(windLen/2);
increment=ceil((1-overlapPct)*windLen);
windMids=round(windLen/2):increment:numel(SUT)-round(windLen/2)-1; % center of each stft window
inds2d=oneWindInd+windMids';
SUT2D=SUT(inds2d').*hann(windLen);
zp=zeros(size(SUT2D));
SUT2D=[zp;zp;SUT2D;zp;zp];
dt=mean(diff(t));
out=fftshift(fft(ifftshift(SUT2D)))*dt;
tout=t(windMids);
fout=linspace(-1/(2*dt),1/(2*dt),windLen*5);


end