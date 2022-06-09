function [phase]=genDiscPhase(lent,level,m,dt,phaseVals,offset)
%level is the time/freq length of each phase level. Offset sets the phase to have the middle of a phaselevel at the center of lent
% 
% 
% if length(phaseVals)~=m
%    'error! not enough phase vals!'
% end

m=numel(phaseVals);

npLevel=round(level/(dt));%Number of points per unit. Unit~Tr contains m pulses, contained within period T1*m. level_t is Tr/m
npUnit=npLevel*m;

nReps=ceil(lent/npUnit);

% unit=zeros(npLevel,m);%Square matrix m by npLevel_t, with all phase values
% for i=1:numel(unit(1,:))
%     unit(:,i)=phaseVals(i);%GaussVals created at begining, from s & m.
%     
% end

unit=repelem(phaseVals,npLevel); %repeat each phase val 



phaseTemporal=repmat(unit,1,nReps); 
phase= phaseTemporal(1:lent);
if offset==1
    
   cent=round(lent/2);
   left=mod(cent,npLevel);
   
   shift=left-round(npLevel/2);
   phase = circshift(phase,shift);
   
end
end
