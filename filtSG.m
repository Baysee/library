function [sigFilt,sig_f,wind,sigFilt_f]=filtSG(sig,nBW,m,displayPlot)
%,varargin)
% Varargin should

lent=numel(sig);
inds=1:lent;


sig_f=nfft(sig,1);

wind=zeros(1,lent);
wind=superGauss(0,nBW,m,inds,round(lent/2));

sigFilt_f=sig_f.*wind;
sigFilt=nifft(sigFilt_f,1);

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