function [xNew,y,shift]=alignXcorr(x,y)


nthpeak=1;
% crossCorr=xcorr(x,y);
if size(x)~=size(y)
    y=y.';
    if size(x)~=size(y)
        error('must have same number of elements')
    end
end
[~,crossCorr]=circCrossCorr(x,y,0);
% [~,locs]=findpeaks(crossCorr,'SortStr','descend','Npeaks',nthpeak);%Find two highest peaks
% shift=-(locs(end)-lent)
[~,locs]=max(crossCorr);
lent=numel(x);
shift=round(locs-lent/2);
xNew=circshift(x,shift);

% figure;plot(xNew); yyaxis right;plot(y); xlim([1.5e4 1.57e4])

