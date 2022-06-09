function [sigFilt,sig_f,wind,sigFilt_f]=filtSG_tf(sig,t,f,nBW,m,displayPlot)
% filtSG_tf(sig,t,f,nBW,m,displayPlot)


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
    plot(f,real(sig_f));
    hold on
    plot(f,real(sigFilt_f));
    yyaxis right
    plot(f,wind);
    subplot(2,1,1)
    plot(t,real(sig))
    hold on
    plot(t,imag(sig))
    plot(t,abs(sig))
    plot(t,real(sigFilt))
    plot(t,imag(sigFilt))
    plot(t,abs(sigFilt))
    legend('real in','imag in','abs in','real out','imag out', 'abs out')
end

end