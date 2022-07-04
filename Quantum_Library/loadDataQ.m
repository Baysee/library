function data=loadDataQ(runDir,filename,lastIndIni,varargin)
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
                    
                    % get a filtered version of the hist to adjust for
                    % temporal drifts
                    indsHist=1:lastIndIni;%19414; 
                    filtBW=2^9; m=8;
                    dat1=data(data(:,1)==1,:);
                    refHist=histcounts(dat1(:,3),0:2^15);
                    refHist=refHist(indsHist);
                    filtRefHist=real(filtSG(refHist,filtBW,m,0));
                else
                    % Adjust for temporal shift
                    dat1_new=data(data(:,1)==1,:);
                    testHist=histcounts(dat1_new(:,3),0:2^15);
                    testHist=testHist(indsHist);
                    filtTestHist=real(filtSG(testHist,filtBW,m,0));

                    [~,~,shift]=alignXcorr( filtTestHist,filtRefHist);
%                      figure;plot(filtRefHist)
%                        hold on
%                     plot(filtTestHist); 
%                     plot(filtTestHistAligned);
                    data(:,3)=mod(data(:,3)+shift,lastIndIni);

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