function cleanedup_angle=clean_angle(sig)
%this function samples the angle to points with a non-zero abs value.
%The angle is then restricted to the 0 to 2*pi domain
bare_angle=angle(sig);
%Sampled_sig returns a signal with a value of one whereever sig>1
tol=0.2;
sampled_sig=abs(sig)>tol*max(abs(sig));
%select only the angle where sig>tol*maxVal
sampled_angle=sampled_sig.*bare_angle;

% if min(sampled_angle<0)
%     sampled_angle=sampled_angle+2*pi;%make all values positive. 
% end
cleanedup_angle=mod(sampled_angle,(2*pi));
