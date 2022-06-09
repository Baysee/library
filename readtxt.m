function [xdat,ydat]=readtxt(filename)%,keepOffset,keepdB)


fid=fopen(filename);
%ydat=10^(ydat/10);
    datacell=textscan(fid,'%f%f','HeaderLines',3);
    xdat=datacell{1};
    ydat=datacell{2};
%     ydat=10.^(ydat/10);


%     if ~exist('keepdB','var')fi
%         ydat=10.^(ydat/10);
%     end
%     
%     if ~exist('keepOffset','var')
%         avgRegion=2000;
%     ydat=ydat-mean([ydat(36924:38000)]);%mean([ydat(1:avgRegion);ydat(end-avgRegion:end)]);
%     ydat=abs(ydat);
%     else
%     msg='Keeping offset';
%     end


end