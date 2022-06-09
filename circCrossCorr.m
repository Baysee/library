function [y,crossCorr]=circCrossCorr(x1,x2,norm);



len=numel(x1);

if len~=numel(x2)
    error('x1 and x2 must have same length!')
end

y=zeros(size(x1));
crossCorr=y;

for i=1:len
    shifted_x1=circshift(x1,-round(len/2)+i);
    mult=x2.*shifted_x1;
    y=y+mult;
    crossCorr(i)=sum(mult);
    
end
if norm==1
    y=y/max(y)*max(x1);
end

a=1;
