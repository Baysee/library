function [out,offset,scaleFactor]=normalizeSignal(in,nBins,offsetRange,sFRange)

% [counts,centers]=hist(in,nBins);
% 
% fGauss2=fit(centers',counts','gauss2');

% offset=min([fGauss2.b1,fGauss2.b2]); 
% scaleFactor=max([fGauss2.b1,fGauss2.b2])-offset;
% 
% offsetRange=1:2000; sFRange=6000:10000;
offset=mean(in(offsetRange));
scaleFactor=mean(in(sFRange))-offset;


out=(in-offset)/scaleFactor;

figure;plot(out); hold on; plot(offsetRange,out(offsetRange));
plot(sFRange,out(sFRange));
nrmdZero=mean(out(offsetRange))
nrmdOne=mean(out(sFRange))
% yyaxis right
% plot(out)
end