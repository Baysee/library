function y=circAutoCorr(x);

len=numel(x);

y=zeros(len,1);

xi=x;

for i=1:len
    
    y(i)=sum(xi.*circshift(x,i));

end

