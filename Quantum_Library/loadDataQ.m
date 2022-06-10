function data=loadDataQ(runDir,filename,varargin)
% Read quantum data when there are many files.
% Make varargin exist to force the function to read only one file

load([runDir filename '.mat']);% dirDat basefilename  '_' num2str(fileNums-1) '.mat'])
data(data(:,1)==0,:)=[];
lastGarbageData=find(data(1:3e3,1)>4,1,'last');
data(1:(lastGarbageData+1),:)=[];

if numel(varargin)>0
    return
end

nMeas=numel(dir([runDir filename(1:end-3) ,'*']));

if nMeas>1
    clockGap=2; % Put this amount of clock triggers between two subsequent measurement. Large enough to not account for a coincidence (i.e., bigger than trig delay) but not big enough to significantly alter total integration time
    for iFn_fileMeas=1:nMeas
        
if iFn_fileMeas==1
    dataAll=data;
else
    data(:,2)=data(:,2)-min(data(:,2))+dataAll(end,2)+clockGap;
    dataAll=[dataAll;data];
end
    end
end

data=dataAll;% clear dataAll

end