function data=loadDataQ(runDir,filename,varargin)
% Read quantum data when there are many files.
% Make varargin exist and equal 1 to force the function to read only one file


if numel(varargin)>0
    if varargin{1}
        load([runDir filename '.mat']);% dirDat basefilename  '_' num2str(fileNums-1) '.mat'])
        data(data(:,1)==0,:)=[];
        lastGarbageData=find(data(1:3e3,1)>4,1,'last');
        data(1:(lastGarbageData+1),:)=[];
        return
    end
end

nMeas=numel(dir([runDir filename(1:end-3) ,'*_0.mat']));
dataAll=[];
% if nMeas>1
    clockGap=2; % Put this amount of clock triggers between two subsequent measurement. Large enough to not account for a coincidence (i.e., bigger than trig delay) but not big enough to significantly alter total integration time
    for iFn_fileMeas=1:nMeas
        
        continueMeasuring=1;
        fileIndex=0;
        try
            while continueMeasuring
                
                
                
                % load data from other files
                load([runDir filename(1:end-3)  num2str(iFn_fileMeas) '_' num2str(fileIndex) '.mat']);% dirDat basefilename  '_' num2str(fileNums-1) '.mat'])
                data(data(:,1)==0,:)=[];
                lastGarbageData=find(data(1:3e3,1)>4,1,'last');
                data(1:(lastGarbageData+1),:)=[];
                
                if iFn_fileMeas==1
                    dataAll=data;
                else
                    data(:,2)=data(:,2)-min(data(:,2))+dataAll(end,2)+clockGap;
                    dataAll=[dataAll;data];
                end
                
                
                fileIndex=fileIndex+1;
                
                if numel(dir([runDir filename(1:end-3) num2str(iFn_fileMeas) '_' num2str(fileIndex) '.mat']))==0
                    continueMeasuring=0;
                end
            end
        catch ME
            % Send e-mail if a crash arises\
            mainError=string(ME.message);
            mestack=string(struct2cell(ME.stack));
            msgToSend={mainError{1},mestack{1},mestack{2},mestack{3}}
            ['Tried to read unexisting file!;'; mainError]
            continue
        end
    end
% end

data=dataAll;% clear dataAll

end