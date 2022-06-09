function [output]=filterSig(m,t0,sig,plotit)

inds=1:numel(sig); half=round(numel(sig)/2);

Sig_f=nrmd_fft(sig,1,1);
filt=superGauss(0,t0,m,inds,half);

sigSize=size(sig);
filtSize=size(filt);

if sigSize(1)~=filtSize
    filt=filt'
end
output=nrmd_ifft(Sig_f.*filt,1,1);
if plotit
    
    figure;
    subplot(2,1,1)
    plot(abs(Sig_f));
    yyaxis right
    plot(filt)
    subplot(2,1,2)
    plot(abs(sig))
    hold on
    plot(abs(output))
    plot(real(output))
    legend('input','abs','real')
end


end


