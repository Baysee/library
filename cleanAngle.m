function [CA]=cleanAngle(y,varargin)

if numel(varargin)==0
    thresh=0.1;
else
    thresh=varargin{1};
end

CA=angle(y);
tooSmall=abs(y)<(thresh*(max(abs(y))-mean(abs(y))));
CA(tooSmall)=0;
CA=unwrap(CA);
