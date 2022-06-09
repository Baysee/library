function [indsOut,pks]=getPeaks( mtest , T_UT )
% mtest is a 2D matrix, with the time and y info. T_UT is the expected separation of the peaks, in units of the time axis in mtest
% basic mtest info
mtestx = mtest(:,1); mtesty = mtest(:,2);
len = numel(mtesty);
dt = mean(diff(mtestx(1:500)));
minDist = round(T_UT/dt*0.99);                  %distance between output peaks

%find peaks by cross correlation
[cc,~] = xcorr(mtesty,(mtesty));
ccL = numel(cc); ccLcenter=round(ccL/2);
cc = cc(ccLcenter:end);% lags=lags(ccLcenter:end);


[peaks,inds]=findpeaks(cc,'MinPeakDistance',minDist);

allinds=zeros(1,len);
allinds(inds)=1;
range=2*minDist;

if len>1.2e6
    topI=1e6;
    mty=mtesty(1:topI);
    ai=allinds(1:topI);
    
    matAI=zeros(2*range,numel(mty));
    sl=-range:range;
%     matAI(range,:)=ai;
    for i=sl
         matAI(i+range+1,:)= circshift(ai,i);
%         matAI(range-i,:)= [ai(i:end), ai(1:i-1)];
%         matAI(i+range,:)= [ai(len-i+1:len), ai(1:len-i)];
    end

    AlignmentVec=matAI*mty;
    [bb,AlInd]=max(AlignmentVec);
 [cc,AlIndInv]=max(-AlignmentVec); % may need to take inverse if
% pulses are inverted!
if cc>bb
    AlInd=AlIndInv;
    warning('Negative pulse was chosen in the function getPeaks!')
end
    allinds=logical(circshift(allinds,sl(AlInd)));
    
else
    
    
    
matAllInds=zeros(2*range,numel(mtesty));
sl=-range:range;
matAllInds(range,:)=allinds;
for i=1:range-1
%     matAllInds(i+range+1,:)= circshift(allinds,i);
    matAllInds(range-i,:)= [allinds(i:end), allinds(1:i-1)];
    matAllInds(i+range,:)= [allinds(len-i+1:len), allinds(1:len-i)];
end

AlignmentVec=matAllInds*mtesty;
[~,AlInd]=max(AlignmentVec);
allinds=logical(matAllInds(AlInd,:));
end



indsOut=1:len; indsOut=indsOut(allinds);

r2=2; inds2test=[-2:2]';
it=inds2test+indsOut;
toobig=it>numel(mtesty); it(toobig)=numel(mtesty);
toosmall=it<1; it(toosmall)=1;
allpks=mtesty(it);
[~,bestpk]=max(allpks);

convBack=[0:2*r2+1:numel(it)-2*r2+1];
indsOut=it(bestpk+convBack);
pks=mtesty(indsOut);
% 
% figure;plot(mtestx,mtesty)
% hold on
% plot(mtestx(indsOut),mtesty(indsOut),'*')
% figure;plot(lags,cc)
% hold on
% plot(inds,peaks,'*')
end