function [freqCent,powerf_linNorm,lambda,power]=getOSAdatBroken(fn,varargin)

[lambda,power]=loadCellData(fn);%[generalFolder, apexFolder, fnFolder, osabasefn]);
% if lambda(end)<lambda(1)
%     lambda=flip(lambda);
%     power=flip(power);
% end

c=299792458;
freqOSA_G=c./lambda;
freqOSA_G=freqOSA_G-mean(freqOSA_G);
% freqOSA_G=fliplr(freqOSA_G);
power_f=flip(power);
powerf_lin=10.^(power_f/10)/max(10.^(power_f/10));
if numel(varargin)
    filtBW=varargin{1};
    plotFiltResults=1;
    superGaussFilt_m=1;
powerf_linNorm=real(filtSG(powerf_lin,filtBW,superGaussFilt_m,plotFiltResults));
powerf_linNorm=powerf_linNorm./max(powerf_linNorm);
else
    powerf_linNorm=powerf_lin./max(powerf_lin);
end
centerFreq=freqOSA_G*powerf_linNorm'/(sum(powerf_linNorm));
freqCent=freqOSA_G-centerFreq;

