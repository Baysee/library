function [sigFiltOut,sig_f,wind,sigFilt_f]=filtSG2D(sigIn,nBW,m,displayPlot)
%,varargin)
% Varargin should

sSig=size(sigIn);
[lent,dim]=max(sSig);
inds=1:lent;

if dim==2
    sig=sigIn.';
else
    sig=sigIn;
end

sig_f=fftshift(fft(ifftshift(sig)));

wind=superGauss(0,nBW,m,inds,round(lent/2))';

sigFilt_f=sig_f.*wind;
sigFilt=fftshift(ifft(ifftshift(sigFilt_f)));

if dim==2
    sigFiltOut=sigFilt.';
else
    sigFiltOut=sigFilt;
end

if displayPlot==1
    figure;
    subplot(2,1,2)
    plot(real(sig_f));
    hold on
    plot(real(sigFilt_f));
    yyaxis right
    plot(wind);
    subplot(2,1,1)
    plot(real(sig))
    hold on
    plot(abs(sig))
    plot(real(sigFilt))
    plot(abs(sigFilt))
    legend('Real input','abs input','real output','abs output')
end


   



end