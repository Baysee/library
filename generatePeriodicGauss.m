function y=generatePeriodicGauss(x,Tr,fwhm,center,m)


sigma=fwhm/(2*sqrt(2*log(2)));

dx=mean(diff(x));
nSingle=round(Tr/dx);
nReps=ceil(numel(x)/nSingle);
single=singleGauss(sigma,0,linspace(-nSingle/2*dx,nSingle/2*dx,nSingle),0);
% single=superGauss(0,sigma,m,linspace(-nSingle/2*dx,nSingle/2*dx,nSingle),0);
% single=sinc(1/sigma*linspace(-nSingle/2*dx,nSingle/2*dx,nSingle));%superGauss(0,sigma,m,linspace(-nSingle/2*dx,nSingle/2*dx,nSingle),0);
y=repmat(single,1,nReps);
y=y(1:numel(x));

if center==1
    lent=numel(x)
   cent=round(lent/2);
   left=mod(cent,nSingle);
   
   shift=left-round(nSingle/2);
   y = circshift(y,shift);
   
end

end