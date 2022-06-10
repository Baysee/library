% get Heralding efficiency, mean photon number, brightness, heralded and unheralded g2
% This program cycles through the different TAI settings (pm on and off),
% and processes the data through various ways (taking only the peaks,
% taking the data between the peaks, etc.). This is done for the different
% voltages, corresponding to different pump powers. Three different noise
% settings are available: no injected noise, 2e4 counts/s of noise, 4e4
% counts/s. 


for FullRuns=2:3
% paths with basic scripts used to process data
addpath('C:\Users\Lord Photon\Documents\MATLAB\Instrumentation\Lab Interfacing\HydraHarp_DataTreatment')
addpath( 'C:\Users\Lord Photon\Documents\MATLAB\library' )

% Directory with all the data
dirDat='C:\Users\Lord Photon\Documents\MATLAB\QTAI\highQ\g2Measurements\';

if FullRuns==1
%%% Case 1: no noise
runDir='march25_g2NoNoise\';
basefilename='QTAI_152-166_184-198_HighQ_0p5nsGauss__';
runFileName='1Int_nonoise_diffPM_15minparts-V';
elseif FullRuns==2
%%% Case 2: 2e4 noise
runDir='march28_g2Noise_2e4\';
basefilename='QTAI_152-166_184-198_HighQ_0p5nsGauss__';
runFileName='3hInt_4ch_2e4Noise_15minParts_V';

elseif FullRuns==3
%%% Case 3: 4e4 noise
runDir='march28_g2Noise_4e4\';
basefilename='QTAI_152-166_184-198_HighQ_0p5nsGauss__';
runFileName='105minInt_4ch_4e4Noise_15minParts_V';
end


vs={'0.174' '0.221' '0.245' '0.268' '0.315' '0.41'};        % Voltages used to control the EDFA power
pows=[10.75 19.55 25.32 30 45 60];                           % power measured from tap after ring drop port

tapLossFac=14.13; %5.75;    % pows above is the power measured on the thorlabs monitor,
                            % after the tap pump power within resonance is estimated to be this much higher than pows.
                            % Taking loss into consideration,: tap 3 dB, WDM: 2 dB,
                            % Ring to fiber Facet: 1.5, ring geometry: 6 dB (from opt.
                            % Express paper)

% Details on experimental conditions
nMode=1; repRate=10.396e6; %tInt=30*60; the Rep rate and integration time is now found automatically from the data set.
tq=1600; % Tai peak separation in ps for q=32, disp 600 km. Used to help find peaks
tsRange=280; % Approximative width of a TAI peak.
pulseTrigDelay=0;%% number of trigs between sig and idler (assuming i.e., sig goes through fiber and appears after a different trig
numTaiDSP=6; % Perform three kinds of data analyses, outlined just below.
doOnlyOnePart=0; % If ==1, then take only one data file subpart instead of the entire data set 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main Script %% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for iTAI=1:numTaiDSP 
for iTAI=1:numTaiDSP
%     iTAI=3;
    
    if iTAI==1
        pm=1; % pm means treat the signal by selecting only the values within TAI peaks
        betweenPeaks=0; % Between peaks is to process the left over counts between the TAI peaks instead of the TAI peaks.
        doPreSelection=1;do2DPostSelection=0; % pre0post1=0 means do a preselection; pre0post1=1 means do a postSelection based on coincidences
        pmStr='PMon'; DSP='peaks'; % Strings to be used to naming the output files
    elseif iTAI==2
        pm=1;        
        betweenPeaks=1; %% must run iTAI==1 before any other settings except from
        doPreSelection=1;do2DPostSelection=0;
        pmStr='PMon'; DSP='betweenPeaks';
    elseif iTAI==3
        pm=0;
        betweenPeaks=0;
        doPreSelection=1;do2DPostSelection=0;
        pmStr='PMoff';DSP='allDatPreselect';
    elseif iTAI==4
        pm=1; 
        doPreSelection=0; do2DPostSelection=0;
        pmStr='PMon'; DSP='PMON_noPreSelect';
    elseif iTAI==5
        pm=0; 
        doPreSelection=0;do2DPostSelection=0;
        pmStr='PMoff'; DSP='PMOFF_noPreSelect';
    elseif iTAI==6
        pm=1; 
        doPreSelection=0;do2DPostSelection=1; % Post selection in the 2D arrivals
        pmStr='PMon'; DSP='PMON_2DPostSelect';
    end
    
    
    for vi=1:numel(vs)
%     for vi=1:numel(vs)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Must be done in two parts. Data was taken in ~15 minutes
        %%% interval, switching between pm on and pm off to minimize
        %%% impacts of potential drifts in the system. Secondly, if the
        %%% .ptu file is too big (>500 Mb), then the mat file is split into
        %%% multiple files
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Determine number of parts (each long measurement is split into multiple
        % parts)
        vStr=char(vs{vi}); % voltage of EDFA in char
        fullBaseName=[ dirDat runDir basefilename runFileName vStr '-' pmStr '-'];
        allParts=dir([fullBaseName '*']); % Find number of parts
        b=struct2cell(allParts);
        fnames_Exp=b(1,:); %% All files that start with the right name
        partCharPos=numel(fullBaseName)-(numel(dirDat)+numel(runDir));
        for fi=1:numel(fnames_Exp)
            fiName=fnames_Exp{fi};
            partNum(fi)=fiName(partCharPos+1);
        end
        nParts=numel(unique(partNum));
        
        % Go through each part (separate measurements I did), as well as 
        %each"subpart", since each ptu file is split into many files when the size exceeds about 500 Mb
        datatemp=[]; data=[]; 
        for ii=1:nParts
            
            % Big files are split into many files with ending _1, _2, ... check if such
            % files exist, if so, append data.
            allmatParts=dir([fullBaseName   num2str(ii) '_*.mat']); % Find number of parts
            b=struct2cell(allmatParts);
            fnames_mat=b(1,:); %% All files that start with the right name
            dirFnames=b(2,:);
            
            for fileNums=1:numel(fnames_mat)
                datatemp=data;
                load([dirFnames{1} '\' fnames_mat{fileNums}]);% dirDat basefilename  '_' num2str(fileNums-1) '.mat'])
                data(data(:,1)==0,:)=[];
                lastGarbageData=find(data(1:3e3,1)>4,1,'last');
                data(1:(lastGarbageData+1),:)=[];
                datatemp=[datatemp;data];
            end
            

            if doOnlyOnePart==1
            % if doOnlyOnePart==1, then get out of nParts loop and proceed with analysis with only one data file
                break; end
            
            if ii~=1
                data(:,2)=data(:,2)-data(1,2)+datatemp(end,2); % Make trigger counts continuous (when different traces are taken experimentally, the trigger count is not continuous)
            end
            
            datatemp=[datatemp;data];
            
        end
        
        data=datatemp; clear datatemp;
        
        %%% This section can be used to reassign certain channels
        % % data(data(:,1)==1,:)=[];
        %  data(data(:,1)==2,1)=3;
        % data(data(:,1)==3,1)=2;
        % data(data(:,1)==4,1)=2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Basic Data Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% This section does basic processing namely assigning data per
        %%% channel, plotting basic histograms and finding the period of
        %%% the SUT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        indsDat=1:numel(data(:,1));                             % Separate data set to each detector. keep track of indeces
        data1=data(data(:,1)==1,:);
        indsDat1=indsDat(data(:,1)==1);
        data2=data(data(:,1)==2,:);
        indsDat2=indsDat(data(:,1)==2);                            % Separate data set to each detector. keep track of indeces
        data3=data(data(:,1)==3,:);
        indsDat3=indsDat(data(:,1)==3);
        data4=data(data(:,1)==4,:);
        indsDat4=indsDat(data(:,1)==4);
        
        [hist1,~,hist1_inds]=histcounts(data1(:,3),0:2^15);     % Sort into histogram. keep track of which index gets assigned to which bin
        [hist2,~,hist2_inds]=histcounts(data2(:,3),0:2^15);
        [hist3,~,hist3_inds]=histcounts(data3(:,3),0:2^15);     % Sort into histogram. keep track of which index gets assigned to which bin
        [hist4,~,hist4_inds]=histcounts(data4(:,3),0:2^15);
        hists=[hist1;hist2;hist3;hist4];
        
        res=4; % HHresolution. Ok to hard code (ps)
        t = (1:2^15)*res;
        lastInd=find(hist1>5,1,'last'); % Last index of the histogram with meaninfull information
        trigPeriod=t(lastInd)*1e-12; % Correponds to period of pulsed signal
        totalIntegrationTime=(data(end,2)-data(1,2))*trigPeriod;
        Np=data(end,2)-data(1,2); % Number of pulses is simply the total number of trigger clocks

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Pre-selection Data Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% If preselection is activated, then only the counts occuring
        %%% within appropriate windows will be kept. Note that here, in
        %%% order to have similar ranges for the PM off case and PM on
        %%% case, the pm off case uses the range from the pm on. If it's
        %%% empty, then it will throw a warning and find a range
        %%% automatically
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        histPks=zeros(4,1);
        
        %%% Tried using filtering... Does not seem useful
        %%% anymore. Keeping here in case I find a use for it.
        %hists_FT=nfft(hists(:,1:lastInd)',1);
        %wind=superGauss(0,3e2,10,1:lastInd,round(lastInd/2));
        %hists_filt=nifft(hists_FT.*wind',1)';
        %hists_filt(:,1:100)=0; hists_filt(:,(lastInd-100):end)=0;
        % figure;plot(abs(hists_filt(1,:)))
        
        if doPreSelection % If doPreSelection==0, then skip this part and use raw data
            
            if pm
                if betweenPeaks==0
                    %% preselection for TAI data if between peaks=0 (so only data within peaks)
                    [taiTW,time_win]= findTAIpks(hists,tq,res,tsRange);
                    taiTW_vi{vi}=taiTW; % Save the indeces of the TAI peaks for this voltage run, containing the indeces for the 4 channels
                else %% if between peaks ==1
                    taiTW=taiTW_vi{vi};
                     for i=1:4                      
                    taiWinds=taiTW{i};
                    %                 if isempty(taiWinds)
                    %                     continue; end
                    fullRange=taiWinds(1):taiWinds(end);
                    fullRange(ismember(fullRange,taiWinds))=[];
                    time_win{i}=fullRange;
                     end
                end
                
            else %% Else if pm=0
                taiTW=taiTW_vi{vi};
                for i=1:4
                    % If available, take data from TAI, with the full range
                    % (peaks and in between peaks)
                    taiWinds=taiTW{i};
                    if isempty(taiWinds) %% This should not happen -- leave here now just to help analyze data temporarily (not for "final" analysis")
                        warning('Need to find data for pm off -- No TAI data available!')
                        noiseFloor(i)=mean(abs(hists_filt(i,200:400)));
                        [histPks(i),histPksLoc(i)]=max(hists_filt(i,:));
                        % Take all points until count is equal to tolNoise*noiseFloor
                        peakMinHeightRatio=0.05;
                        rangeL=find(abs(hists_filt(i,1:histPksLoc(i)))>peakMinHeightRatio*(histPks(i)-noiseFloor(i))+noiseFloor(i),1);
                        rangeR=find(abs(hists_filt(i,histPksLoc(i):end))>peakMinHeightRatio*(histPks(i)-noiseFloor(i))+noiseFloor(i),1,'last');
                        time_win(i)={rangeL:(histPksLoc(i)+rangeR)};
                    else
                        fullRange=taiWinds(1):taiWinds(end);
                        time_win{i}=fullRange;
                    end
                end
            end
            
            
            %% Select points from pre-selection
            [inds1,data1_tm,inds1_tm,hist1_tm]=preSelect(time_win{1},hist1_inds,data1,indsDat1);
            [inds2,data2_tm,inds2_tm,hist2_tm]=preSelect(time_win{2},hist2_inds,data2,indsDat2);
            [inds3,data3_tm,inds3_tm,hist3_tm]=preSelect(time_win{3},hist3_inds,data3,indsDat3);
            [inds4,data4_tm,inds4_tm,hist4_tm]=preSelect(time_win{4},hist4_inds,data4,indsDat4);
            
            [~,order]=sort([inds1_tm,inds2_tm,inds3_tm,inds4_tm]);
            dataSelected=[data1_tm;data2_tm;data3_tm;data4_tm]; dataSelected=dataSelected(order,:);   % Put it back into the order at which it was detected for the "trigger selection"
            
        else
            % If no preselection is done, simply assign all data to the
            % place holders. 
            dataSelected=data;
            hist1_tm=hist1;hist2_tm=hist2;hist3_tm=hist3;hist4_tm=hist4;
        end
        
        
        %% To ease memory load, clear data and data1
        clear data datatemp;
        clear data1 data2 data3 data4
        %% Plot selected regions
        %
        % figure
        % plot(hist1);hold on; plot(hist1_tm)
        % plot(hist2);hold on; plot(hist2_tm)
        % plot(hist3);hold on; plot(hist3_tm)
        % plot(hist4);hold on; plot(hist4_tm)
%         figure
%         plot(t/1e3,hist1);hold on; plot(t/1e3,hist1_tm)
%         plot(t/1e3,hist2);hold on; plot(t/1e3,hist2_tm)
%         plot(t/1e3,hist3);hold on; plot(t/1e3,hist3_tm)
%         plot(t/1e3,hist4);hold on; plot(t/1e3,hist4_tm)
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Get Coincidences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Find double and triple coincidences %%%%%%%%%%%%%%%%%%%%%%%%%
        % function is as follows:
        % [datadouble12_1,datadouble12_2]=findDoubleCoinc(pulseTrigDelay,1,2,dataSelected,plotHist,plotNumEvents);
        % maps to:
        % [outputDataCoincidence_fromChannel_c1,outputDataCoincidence_fromChannel_c2]=...
        % findDoubleCoinc(Number of pulse trig delays that are recognized as a coincidence,...
        % first channel label for c1, second channel c2 ,plotHist=1 to plot histogram of trigg delays
        % ,plotNumEvents=1 to plot the number of coincidences found at each
        % shift); Similarly for triple coincidence function.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        plotHist=0;plotNumEvents=0;
        
        %%% "Self coincidences" (from both outputs of a single BS).
        c1=1;c2=2;
        dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis
        [datadouble12_1,datadouble12_2]=findDoubleCoincSelectedDat(pulseTrigDelay,dataCoinc,plotHist,plotNumEvents);
        
        c1=3;c2=4;
         dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis
        [datadouble34_1,datadouble34_2]=findDoubleCoincSelectedDat(pulseTrigDelay,dataCoinc,plotHist,plotNumEvents);
        
        %%% Coincidences between signal and idler
        %%% For each case, also find CAR
        
        taiTW=taiTW_vi{vi}; % Select TAI indeces for this voltage run
        nac=11; coincidenceCount=zeros(1,4); accidentalsCount=zeros(1,4);

        c1=1;c2=3; taiInds1=taiTW{c1}; taiInds2=taiTW{c2};
        dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis
        [coincidenceCount(1),accidentalsCount(1),hist13_psAll,hist13_psAll_iac,datadouble13_1,datadouble13_2]=getCarHistCoinc(dataCoinc,nac,...
        pulseTrigDelay,do2DPostSelection,taiInds1,taiInds2,plotHist,plotNumEvents);
    
        c1=1;c2=4; taiInds1=taiTW{c1}; taiInds2=taiTW{c2};     
        dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis
        [coincidenceCount(2),accidentalsCount(2),hist14_psAll,hist14_psAll_iac,datadouble14_1,datadouble14_2]=getCarHistCoinc(dataCoinc,nac,...
        pulseTrigDelay,do2DPostSelection,taiInds1,taiInds2,plotHist,plotNumEvents);
    
        c1=2;c2=3; taiInds1=taiTW{c1}; taiInds2=taiTW{c2};
        dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis
        [coincidenceCount(3),accidentalsCount(3),hist23_psAll,hist23_psAll_iac,datadouble23_1,datadouble23_2]=getCarHistCoinc(dataCoinc,nac,...
        pulseTrigDelay,do2DPostSelection,taiInds1,taiInds2,plotHist,plotNumEvents);
    
        c1=2;c2=4; taiInds1=taiTW{c1}; taiInds2=taiTW{c2};   
        dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis
        [coincidenceCount(4),accidentalsCount(4),hist24_psAll,hist24_psAll_iac,datadouble24_1,datadouble24_2]=getCarHistCoinc(dataCoinc,nac,...
        pulseTrigDelay,do2DPostSelection,taiInds1,taiInds2,plotHist,plotNumEvents);

    %%% 2D post selection is done within coincidence function instead
%         if do2DPostSelection
%             taiTW=taiTW_vi{vi}; % Get indeces of TAI peaks for this voltage setting
%             c1=1;c2=3;
%             dd1=datadouble13_1; dd2=datadouble13_2;
%             taiInds1=taiTW{c1}; taiInds2=taiTW{c2};
%             postSelectedInd=ismember(dd1(:,3),taiInds1) | ismember(dd2(:,3),taiInds2) ;
%             dd1Post=dd1(postSelectedInd,:);            dd2Post=dd2(postSelectedInd,:);
%         end
%     
    CC=sum(coincidenceCount); AC=sum(accidentalsCount); CAR=(CC-AC)/AC; dCAR=CAR*sqrt(1/CC+1/AC);


%         c1=1;c2=4;
%         [datadouble14_1,datadouble14_2]=findDoubleCoinc(pulseTrigDelay,1,4,dataSelected,plotHist,plotNumEvents);
%         c1=2;c2=3;
%         [datadouble23_1,datadouble23_2]=findDoubleCoinc(pulseTrigDelay,2,3,dataSelected,plotHist,plotNumEvents);
%         c1=2;c2=4;
%         [datadouble24_1,datadouble24_2]=findDoubleCoinc(pulseTrigDelay,2,4,dataSelected,plotHist,plotNumEvents);
%         
        
        %% Plot correlated counts on t1-t2 axis
         timeCoinc = ((1:2^16)-2^15)*res;
         
         %%% Perhaps not much point in getting delays directly with all
         %%% channels?
%          photDet1=[datadouble13_1;datadouble14_1;datadouble23_1;datadouble24_1]; % Channel 1 and 2
%          photDet2=[datadouble13_2;datadouble14_2;datadouble23_2;datadouble24_2]; % Channel 3 and 4 
%         [hist1_psAll,~,hist1_psAll_inds]=histcounts(photDet1(:,3),0:2^15);     % Sort into histogram. keep track of which index gets assigned to which bin
%         [hist2_psAll,~,hist2_psAll_inds]=histcounts(photDet2(:,3),0:2^15); 
%        [hist1234_psAll,~,hist1234_psAll_inds]=histcounts(photDet1(:,3)-photDet2(:,3)+2^15,0:2^16);


        % get histogram of time delays
%         [hist13_psAll,~,hist13_psAll_inds]=histcounts(datadouble13_1(:,3)-datadouble13_2(:,3)+2^15,0:2^16);
%         [hist14_psAll,~,hist14_psAll_inds]=histcounts(datadouble14_1(:,3)-datadouble14_2(:,3)+2^15,0:2^16);
%         [hist23_psAll,~,hist23_psAll_inds]=histcounts(datadouble23_1(:,3)-datadouble23_2(:,3)+2^15,0:2^16);
%         [hist24_psAll,~,hist24_psAll_inds]=histcounts(datadouble24_1(:,3)-datadouble24_2(:,3)+2^15,0:2^16);

%         figure;subplot(2,1,1)
%         plot(timeCoinc,hist13_psAll)
%         hold on; plot(timeCoinc,hist14_psAll)
%         plot(timeCoinc,hist23_psAll); plot(timeCoinc,hist24_psAll)
%         xlabel('time (ps)');ylabel('G_2(\tau) ');legend('ch1-ch3','ch1-ch4','ch2-ch3','ch2-ch4')
%         title(['CAR: ' num2str(CAR)])
%         subplot(2,1,2)
%         plot(timeCoinc,hist13_psAll_iac)
%         hold on; plot(timeCoinc,hist14_psAll_iac)
%         plot(timeCoinc,hist23_psAll_iac); plot(timeCoinc,hist24_psAll_iac)
%         xlabel('time (ps)');ylabel('G_2(\tau) (accidentals)');legend('ch1-ch3','ch1-ch4','ch2-ch3','ch2-ch4')
%         allCoinc=numel(datadouble13_1(:,3))+numel(datadouble14_1(:,3))+numel(datadouble23_1(:,3))+numel(datadouble24_1(:,3))
                
        %%% Get 2D representation
%         datadouble13_2D=get2Drep(datadouble13_1(:,3),datadouble13_2(:,3)); plotReduced2Drep(datadouble13_2D,hist1,hist3)
%         datadouble14_2D=get2Drep(datadouble14_1(:,3),datadouble14_2(:,3));
%         datadouble23_2D=get2Drep(datadouble23_1(:,3),datadouble23_2(:,3));
%         datadouble24_2D=get2Drep(datadouble24_1(:,3),datadouble24_2(:,3));
%         

        %% Plot hists of coincidences
        
        data1_cc=[datadouble13_1;datadouble14_1];
        data2_cc=[datadouble23_1;datadouble24_1];
        data3_cc=[datadouble13_2;datadouble23_2];
        data4_cc=[datadouble14_2;datadouble24_2];
        
        [hist1_cc,~,hist1_inds_cc]=histcounts(data1_cc(:,3),0:2^15);     % Sort into histogram. keep track of which index gets assigned to which bin
        [hist2_cc,~,hist2_inds_cc]=histcounts(data2_cc(:,3),0:2^15);
        [hist3_cc,~,hist3_inds_cc]=histcounts(data3_cc(:,3),0:2^15);     % Sort into histogram. keep track of which index gets assigned to which bin
        [hist4_cc,~,hist4_inds_cc]=histcounts(data4_cc(:,3),0:2^15);
        
        
        %% Triple Coincidences
        c1=1;c2=2;c3=3;
        [dataTrip123_1,dataTrip123_2,dataTrip123_3]=findTripleCoinc(pulseTrigDelay,c1,c2,c3,dataSelected,plotNumEvents);
        c1=1;c2=2;c3=4;
        [dataTrip124_1,dataTrip124_2,dataTrip124_3]=findTripleCoinc(pulseTrigDelay,c1,c2,c3,dataSelected,plotNumEvents);
        c1=1;c2=3;c3=4;
        [dataTrip134_1,dataTrip134_2,dataTrip134_3]=findTripleCoinc(pulseTrigDelay,c1,c2,c3,dataSelected,plotNumEvents);
        c1=2;c2=3;c3=4;
        [dataTrip234_1,dataTrip234_2,dataTrip234_3]=findTripleCoinc(pulseTrigDelay,c1,c2,c3,dataSelected,plotNumEvents);
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Get all counts and calculate metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        
        
        n1=sum(hist1_tm) % Per channel, total number of counts after preselection
        n2=sum(hist2_tm)
        n3=sum(hist3_tm)
        n4=sum(hist4_tm)
        
        
        n1Raw=sum(hist1); % Per channel, total number of raw counts
        n2Raw=sum(hist2);
        n3Raw=sum(hist3);
        n4Raw=sum(hist4);
        
        c12=numel(datadouble12_1(:,1)) % coincidences between the same beam splitter
        c34=numel(datadouble34_1(:,1))
        
        c13=numel(datadouble13_1(:,1))
        c14=numel(datadouble14_1(:,1))
        c23=numel(datadouble23_1(:,1))
        c24=numel(datadouble24_1(:,1))
        
        c134=numel(dataTrip134_1(:,1))
        c234=numel(dataTrip234_1(:,1))
        c123=numel(dataTrip123_1(:,1))
        c124=numel(dataTrip124_1(:,1))
        
        
        cis=(c13+c14+c23+c24) % Total number of coincidences
        ns=(n1+n2) % Total number of input counts, idler (after preselection, if activated)
        ni=(n3+n4) % Total number of input counts, signal (after preselection, if activated)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Calculate metrics         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% heralded g2
        g2h_i=(c134+c234)*(ns)/((c13+c23)*(c14+c24)),%(c134+c234)*Np/((c13+c23)*(c14+c24))
        dg2h_i=g2h_i*sqrt(1/(c134+c234)+1/(ns)+1/(c13+c23)+1/(c14+c24)),%g2h_i*sqrt(1/(c123+c234)+1/(c13+c23)+1/(c14+c24))
        
        
        g2h_s=(c123+c124)*(ni)/((c13+c14)*(c23+c24))
        dg2h_s=g2h_s*sqrt(1/(c123+c124)+1/(c13+c14)+1/(c23+c24)+1/(ni))
        
        %% unheralded g2
        g2_12=c12*Np/(n1*n2)
        dg2_12=g2_12*sqrt(1/c12+1/n1+1/n2)
        K12=1/(g2_12-1)
        
        g2_34=c34*Np/(n3*n4)
        dg2_34=g2_34*sqrt(1/c34+1/n3+1/n4)
        K34=1/(g2_34-1)
        
        
        
        
        %% heralding eff and mean photon
        
        % g2_12_34=cis*Np/(ni*ni)
        
        etas=cis/ni
        detai=etas*sqrt(1/cis+1/ni)
        etai=cis/ns
        detas=etai*sqrt(1/cis+1/ns)
        
        BrightUnNormed=(n1+n2)*(n3+n4)/(cis*totalIntegrationTime) % Pairs/second per mode
        dBrightUnNormed=BrightUnNormed*sqrt(1/(n1+n2)+1/(n3+n4)+1/cis)
        mpn=(ni)*(ns)/(cis*Np)
        dmpn=mpn*sqrt(1/(n1+n2)+1/(n3+n4)+1/cis)
        
        
        CAR_vi(vi)=CAR; dCAR_vi(vi)=dCAR;
        etas_vi(vi)=etas;detas_vi(vi)=detas;
        etai_vi(vi)=etai;detai_vi(vi)=detai;
        BrightUnNormed_vi(vi)=BrightUnNormed;dBrightUnNormed_vi(vi)=dBrightUnNormed;
        mpn_vi(vi)=mpn;dmpn_vi(vi)=dmpn;
        g2h_i_vi(vi)=g2h_i;dg2h_i_vi(vi)=dg2h_i;
        g2h_s_vi(vi)=g2h_s;dg2h_s_vi(vi)=dg2h_s;
        g2_12_vi(vi)=g2_12;dg2_12_vi(vi)=dg2_12;
        g2_34_vi(vi)=g2_34;dg2_34_vi(vi)=dg2_34;
        
        nCh_vi(vi,:)=[n1,n2,n3,n4];
        coinc_vi(vi,:)=[c13,c14,c23,c24];
        selfCoinc_vi(vi,:)=[c12,c34];
        triples_vi(vi,:)=[c123,c124,c134,c234];
        nRaw_vi(vi,:)=[n1Raw,n2Raw,n3Raw,n4Raw]
        
        
        
        
        %% Save data, metrics, counts.

        
       save([dirDat runDir 'data-' basefilename vStr DSP '.mat'],...'datadouble13_2D','datadouble14_2D','datadouble23_2D','datadouble24_2D',...
            'datadouble13_1','datadouble14_1','datadouble23_1','datadouble24_1',...
            'datadouble13_2','datadouble14_2','datadouble23_2','datadouble24_2',...
            'hist13_psAll','hist14_psAll','hist23_psAll','hist24_psAll',...
            'hist1','hist2','hist3','hist4',...
            'hist1_tm','hist2_tm','hist3_tm','hist4_tm',...
            'timeCoinc','t');
        save([dirDat runDir 'counts-' basefilename vStr DSP '.mat'],'n1','n2','n3','n4',...
            'c13','c14','c23','c24',...
            'c134','c234','c123','c124',...
            'c12','c34','totalIntegrationTime',...
            'cis','ni','ns','Np');
        save([dirDat runDir 'metrics-' basefilename vStr DSP '.mat'],'g2h_i','dg2h_i','g2h_s','dg2h_s',...
            'g2_12','dg2_12','g2_34','dg2_34',...
            'etai','detai','etas','detas','CAR','dCAR',...
            'BrightUnNormed','dBrightUnNormed','mpn','dmpn')
        
        
        h1=figure;
        plot(t/1e3,hist1)
        hold on
        plot(t/1e3,hist2)
        plot(t/1e3,hist3)
        plot(t/1e3,hist4)
        plot(t/1e3,hist1_tm)
        plot(t/1e3,hist2_tm)
        plot(t/1e3,hist3_tm)
        plot(t/1e3,hist4_tm)
        xlabel('Time (ns)');ylabel('Counts')
        set(findall(gcf,'-property','FontSize'),'FontSize',16)
        %     end
        % end
        savefig(h1,[dirDat runDir 'fig-' basefilename vStr pmStr '-' DSP  '.fig']);
        
        
        h2=figure;
        subplot(2,1,1)
        plot(timeCoinc,hist13_psAll)
        hold on; plot(timeCoinc,hist14_psAll)
        plot(timeCoinc,hist23_psAll); plot(timeCoinc,hist24_psAll)
        xlabel('time (ps)');ylabel('G_2(\tau) ');legend('ch1-ch3','ch1-ch4','ch2-ch3','ch2-ch4')
        title(['CAR: ' num2str(CAR)])
        subplot(2,1,2)
        plot(timeCoinc,hist13_psAll_iac)
        hold on; plot(timeCoinc,hist14_psAll_iac)
        plot(timeCoinc,hist23_psAll_iac); plot(timeCoinc,hist24_psAll_iac)
        xlabel('time (ps)');ylabel('G_2(\tau) (accidentals)');legend('ch1-ch3','ch1-ch4','ch2-ch3','ch2-ch4')
%         figure;
%         plot(timeCoinc/1e3,hist12_psAll)
%         xlabel('Time (ns)');ylabel('Counts')
%         set(findall(gcf,'-property','FontSize'),'FontSize',16)
        savefig(h2,[dirDat runDir 'fig-' basefilename vStr pmStr '-' DSP  'CoincDelays.fig']);
        
        close(h1);
        close(h2);
        
    end
    
    
  
    
    
    if iTAI==1
        newFig=1;hetai=[]; hetas=[];hmpn=[]; hBright=[];
        hnRaw=[]; hnCh=[]; hcoinc=[]; hCAR=[];
        saveIt=0;  hg2h_i_vi=[]; hg2h_s_vi=[]; 
        hg2_12_vi=[]; hg2_34_vi=[];
    elseif iTAI==numTaiDSP
        newFig=0; saveIt=1; % On last pass, save figure;
    else
        newFig=0; saveIt=0;
    end
    
    hCAR=plotMetric(pows,CAR_vi,dCAR_vi,newFig,hCAR,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    hetas=plotMetric(pows,etas_vi,detas_vi,newFig,hetas,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    hetai=plotMetric(pows,etai_vi,detai_vi,newFig,hetai,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    hmpn=plotMetric(pows,mpn_vi,dmpn_vi,newFig,hmpn,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    NormalizedBrightness=BrightUnNormed_vi./(pows*tapLossFac);
    dNormalizedBrightness=dBrightUnNormed_vi./(pows*tapLossFac);
    hBright=plotMetric(pows,NormalizedBrightness,dNormalizedBrightness,newFig,hBright,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    
    hnRaw=plotMetric(pows,nRaw_vi,sqrt(nRaw_vi),newFig,hnRaw,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    hnCh=plotMetric(pows,nCh_vi,(sqrt(nCh_vi)),newFig,hnCh,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    hcoinc=plotMetric(pows,coinc_vi,sqrt(coinc_vi),newFig,hcoinc,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    
    hg2h_i_vi=plotMetric(pows,g2h_i_vi,dg2h_i_vi,newFig,hg2h_i_vi,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    hg2h_s_vi=plotMetric(pows,g2h_s_vi,dg2h_s_vi,newFig,hg2h_s_vi,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    hg2_12_vi=plotMetric(pows,g2_12_vi,dg2_12_vi,newFig,hg2_12_vi,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    hg2_34_vi=plotMetric(pows,g2_34_vi,dg2_34_vi,newFig,hg2_34_vi,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
    
    
    nCh_vi(vi,:)=[n1,n2,n3,n4];
    coinc_vi(vi,:)=[c13,c14,c23,c24];
    selfCoinc_vi(vi,:)=[c12,c34];
    triples_vi(vi,:)=[c123,c124,c134,c234];
    nRaw_vi(vi,:)=[n1Raw,n2Raw,n3Raw,n4Raw];
    
    
end





end

function [coincidenceCount,accidentalsCount,hist13_psAll,hist13_psAll_iac,datadouble13_1,datadouble13_2]=getCarHistCoinc(dataCoinc,nac,...
    pulseTrigDelay,postSel,taiInds1,taiInds2,plotHist,plotNumEvents)
%    c1=1;c2=3;        
%         dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis

%         coincidenceCount=zeros(1,4); accidentalsCount=zeros(1,4);

%         nac=11;       
        acc=zeros(1,nac);
        for iac=1:nac
        pulseTrigDelay_iac=iac-1+pulseTrigDelay;
        [datadouble13_1_iac,datadouble13_2_iac]=findDoubleCoincSelectedDat(pulseTrigDelay_iac,dataCoinc,plotHist,plotNumEvents);
        if postSel %% Post selection only if one of two indeces is within the "grid" (one count belongs to a TAI peak)
            postSelectedInd=ismember(datadouble13_1_iac(:,3),taiInds1) | ismember(datadouble13_2_iac(:,3),taiInds2) ;
            datadouble13_1_iac=datadouble13_1_iac(postSelectedInd,:);        
            datadouble13_2_iac=datadouble13_2_iac(postSelectedInd,:);
        end
        
        acc(iac)= numel(datadouble13_1_iac(:,1));
        if iac==1
           [hist13_psAll,~,~]=histcounts(datadouble13_1_iac(:,3)-datadouble13_2_iac(:,3)+2^15,0:2^16);      
            datadouble13_1=datadouble13_1_iac;
            datadouble13_2=datadouble13_2_iac;
        end
        if iac==nac
           [hist13_psAll_iac,~,~]=histcounts(datadouble13_1_iac(:,3)-datadouble13_2_iac(:,3)+2^15,0:2^16);      
        end
        end
        coincidenceCount(1)=acc(1); accidentalsCount(1)=mean(acc(2:end));
%         CAR=(coincidenceCount(1)-accidentalsCount(1))/accidentalsCount(1)

end


function [datadouble2D]=get2Drep(datadouble1,datadouble2)

     datadouble2D=zeros(2^15);
        for i=1:numel(datadouble1)
            event_t1All=datadouble1(i)+1;
            event_t2All=datadouble2(i)+1;

            datadouble2D(event_t1All,event_t2All)=datadouble2D(event_t1All,event_t2All)+1;
        end

end

function plotReduced2Drep(datadouble2D,hist1,hist2)
        range=-1000:999;
        [~,maxIndHist1]=max(hist1); [~,maxIndHist2]=max(hist2);
        datadouble2D_reduced=datadouble2D(maxIndHist1+range,maxIndHist2+range);
        figure;imagesc(datadouble2D_reduced)
        
        
end


function [inds1,data1_tm,inds1_tm,hist1_tm]=preSelect(tm1,hist1_inds,data1,indsDat1)
% tm1 = time_win{1};                                      % select counts within windows
% tm1(tm1<1)=[];
inds1 = ismember(hist1_inds,tm1);                       % Take only the counts there are within the preselections
data1_tm=data1(inds1,:);
inds1_tm=indsDat1(inds1);
[hist1_tm,~]=histcounts(data1_tm(:,3),0:2^15);
end



function [occurence]=findReptitionNumber(y)
% example:
%y=[ 1 2 3 4 4 4 5 5 6]
% occurence = [ 1     1     1     3     2     1]
f=unique(y);
occurence=histc(y,f);
end

function hF=plotMetric(pows,dats,errors,newFig,hF,saveIt,pmState,dirDat,basefilename)

a=inputname(2);
%
% varPlot='BrightUnNormed'
% varPlotErr=['d' varPlot]

if newFig
    hF=figure;
else
    figure(hF)
    hold on
end

nDatToPlot=size(dats,1);
if nDatToPlot>1
    for i=1:size(dats,2)
        if sum((dats(:,i)))==0
            continue
        else
            errorbar(pows,dats(:,i),errors(:,i),'DisplayName',[pmState])
            hold on
        end
    end
else
    errorbar(pows,dats,errors,'DisplayName',[pmState])
end
legend('Show')
xlabel('Power (\mu W)'); ylabel(a)
set(findall(gcf,'-property','FontSize'),'FontSize',16)

if saveIt
    savefig(hF,[dirDat 'MetricsFig-' basefilename a '.fig']);
end
end


function [taiTW,time_win]= findTAIpks(hists,tq,res,tsRange)
nhists=numel(hists(:,1));
for i=1:nhists % cycle through the four channels
    %                 noiseFloor(i)=mean(hists_filt(i,200:400)); % At some
    %                 point, filtering seems useful... not used anymore
    noiseFloor(i)=mean(hists(i,200:400)); % Find noise floor to help find peaks
    peakMinHeightRatio=0.04;
    pkSeparation=round(tq/res*0.5); % Give tolerance to peak separation
    pkWidth=round(tsRange/res);
    pkRange=round(-pkWidth/2):1:round(pkWidth/2);
    [histPks(i),~]=max(hists(i,:));
    [~,locs,w]=findpeaks(hists(i,:),'MinPeakHeight',...
        round(peakMinHeightRatio*(histPks(i)-noiseFloor(i))+noiseFloor(i)),...
        'MinPeakDistance',pkSeparation);
    twi=[];
    for ipw=1:numel(locs)
        twi =[twi,locs(ipw)+(-round(w(ipw)*1.5):round(w(ipw)*1.5))];
    end
    time_win(i)={twi};
    %     time_win(i)={reshape(locs+pkRange',numel(locs)*numel(pkRange),1)};
end
taiTW=time_win;

end

% 
% 
% function [coincidenceCount,accidentalsCount,hist13_psAll,hist13_psAll_iac,datadouble13_1,datadouble13_2]=getCarHistCoinc(dataCoinc,nac,...
%     pulseTrigDelay,plotHist,plotNumEvents)
% %    c1=1;c2=3;        
% %         dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis
% 
% %         coincidenceCount=zeros(1,4); accidentalsCount=zeros(1,4);
% 
% %         nac=11;       
%         acc=zeros(1,nac);
%         for iac=1:nac
%         pulseTrigDelay_iac=iac-1+pulseTrigDelay;
%         [datadouble13_1_iac,datadouble13_2_iac]=findDoubleCoincSelectedDat(pulseTrigDelay_iac,dataCoinc,plotHist,plotNumEvents);
%         acc(iac)= numel(datadouble13_1_iac(:,1));
%         if iac==1
%            [hist13_psAll_iac,~,~]=histcounts(datadouble13_1_iac(:,3)-datadouble13_2_iac(:,3)+2^15,0:2^16);      
%             hist13_psAll=hist13_psAll_iac;
%             datadouble13_1=datadouble13_1_iac;
%             datadouble13_2=datadouble13_2_iac;
%         end
%         end
%         coincidenceCount(1)=acc(1); accidentalsCount(1)=mean(acc(2:end));
% %         CAR=(coincidenceCount(1)-accidentalsCount(1))/accidentalsCount(1)
% 
% end