addpath('/Users/ben/Documents/MATLAB/Instrumentation/Lab Interfacing/HydraHarp_DataTreatment')
addpath( '/Users/ben/Documents/MATLAB/library' )

% Directory with all the data
dirDat='';


runDir='SNSPD_traces/';'sourceBypassed/';'Measurements_25052021-27052021/Port11-0/26052021/';
% basefilename='5-6_sinc_opt_500ps_m25dB_47-148_221-322_4chg2_5min_res4_0';
%  basefilename='power_m8p3dBm_35-39-76_327-331_0';
basefilename='4-7_sinc_opt_500ps_m30dB_142-143_218-219_0';
pulseTrigDelay=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Script %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




load([runDir basefilename '.mat']);% dirDat basefilename  '_' num2str(fileNums-1) '.mat'])
data(data(:,1)==0,:)=[];
lastGarbageData=find(data(1:3e3,1)>4,1,'last');
data(1:(lastGarbageData+1),:)=[];

%%% This section can be used to reassign certain channels
% % data(data(:,1)==1,:)=[];
%  data(data(:,1)==2,1)=3;
% data(data(:,1)==3,1)=2;
% data(data(:,1)==4,1)=2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Data Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This section does basic processing namely assign1ng data per
%%% channel, plotting basic histograms and finding the period of
%%% the SUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


indsDat=1:numel(data(:,1));                             % Separate data set to each detector. keep track of indeces
data1=data(data(:,1)==1,:);
indsDat1=indsDat(data(:,1)==1);
data2=data(data(:,1)==2,:);
indsDat2=indsDat(data(:,1)==2);   
indsDat=1:numel(data(:,1));                             % Separate data set to each detector. keep track of indeces
data3=data(data(:,1)==3,:);
indsDat3=indsDat(data(:,1)==3);
data4=data(data(:,1)==4,:);
indsDat4=indsDat(data(:,1)==4);     % Separate data set to each detector. keep track of indeces


[hist1,~,hist1_inds]=histcounts(data1(:,3),0:2^15);     % Sort into histogram. keep track of which index gets assigned to which bin
[hist2,~,hist2_inds]=histcounts(data2(:,3),0:2^15);
[hist3,~,hist3_inds]=histcounts(data3(:,3),0:2^15);     % Sort into histogram. keep track of which index gets assigned to which bin
[hist4,~,hist4_inds]=histcounts(data4(:,3),0:2^15);
hists=[hist1;hist2];

res=4;256; % HHresolution. Ok to hard code (ps)
t = (1:2^15)*res;
lastInd=find(hist1>5,1,'last'); % Last index of the histogram with mean1nfull information
trigPeriod=t(lastInd)*1e-12; % Correponds to period of pulsed signal
totalIntegrationTime=(data(end,2)-data(1,2))*trigPeriod;
Np=data(end,2)-data(1,2); % Number of pulses is simply the total number of trigger clocks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-selection Data Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

histPks=zeros(4,1);


preSelectionRange_Left=80;  preSelectionRange_Right=300;  
range=-preSelectionRange_Left:preSelectionRange_Right;
[maxVal1,maxInd1]=max(hist1); time_win{1}=maxInd1+range;
[maxVal2,maxInd2]=max(hist2); time_win{2}=maxInd2+range;

[inds1,data1_tm,inds1_tm,hist1_tm]=preSelect(time_win{1},hist1_inds,data1,indsDat1);
[inds2,data2_tm,inds2_tm,hist2_tm]=preSelect(time_win{2},hist2_inds,data2,indsDat2);

[~,order]=sort([inds1_tm,inds2_tm]);%,inds3_tm,inds4_tm]);
dataSelected=[data1_tm;data2_tm];%;data3_tm;data4_tm]; 
dataSelected=dataSelected(order,:);   % Put it back into the order at which it was detected for the "trigger selection"


% a=(dataSelected(:,2)-circshift(dataSelected(:,2),1));
% figure;histogram(a(abs(a)<1e3),1000);

%% Plot selected region2
%
figure
plot(t,hist1);hold on; plot(t,hist1_tm)
plot(t,hist2);hold on; plot(t,hist2_tm)
plot(t,hist3);hold on;% plot(t,hist1_tm)
plot(t,hist4);hold on; %plot(t,hist2_tm)

figure
plot(hist1);hold on; plot(hist1_tm)
plot(hist2);hold on; plot( hist2_tm)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Coincidences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find double and triple coincidences %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% If you want to plot the histogram of the trig delays and the number of
% events per trig delay shift
plotHist=1;plotNumEvents=1;

c1=1;c2=2;

dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis
[datadouble12_1,datadouble12_2,data_photdet1,data_photdet2]=getDoubleCoincUpdated(pulseTrigDelay,dataCoinc,plotHist,plotNumEvents);


%% In order to calculate CAR, need to first find coincidences and remove them from our data set to find the accidentals.
%% %%%%% Plot correlated counts on t1-t2 axis

time12 = ((1:2^16)-2^15)*res;

% "All" in hist12_psAll is before post selection window.
[hist12_psAll,~]=histcounts(datadouble12_1(:,3)-datadouble12_2(:,3)+2^15,0:2^16);
timeCoinc = ((1:2^16)-2^15)*res;
% figure;plot(timeCoinc,hist12_psAll);


[~,maxCorrInd]=max(hist12_psAll);

%% %%%%%% RangeCorrWind is the variable that may needs adjusting! It is related to the time_win variables above, but acts on the correlated counts

rangeCorrWind=150; corrWind=(-rangeCorrWind:rangeCorrWind)+maxCorrInd;


%% now that correlation window is defined, select only those counts.
alltDiffs=datadouble12_1(:,3)-datadouble12_2(:,3)+2^15;
corrCounts=ismember(alltDiffs,corrWind); %find correlated counts that are within the correlated window corrWind

data_photdet1cc=data_photdet1(corrCounts);      % Keep only counts within correlation window
data_photdet2cc=data_photdet2(corrCounts);
datadouble12_1cc=datadouble12_1(corrCounts,:);
datadouble12_2cc=datadouble12_2(corrCounts,:);


[hist12_ps,~]=histcounts(datadouble12_1cc(:,3)-datadouble12_2cc(:,3)+2^15,0:2^16);
figure;plot(timeCoinc,hist12_psAll); hold on; plot(timeCoinc,hist12_ps);


if numel(data_photdet1cc)~=numel(unique(data_photdet1cc)) || numel(data_photdet2cc)~=numel(unique(data_photdet2cc))
 % This should not happen since width of pulses is smaller than dead
 % time... leave rest of code here for later in case this problem arises.
        warning('WARNING: some photon counts appear more than once in the coincidences!!')
        pause
%     [a,b] = histc(phot1_id,unique(phot1_id)); y1 = a(b); y1=y1>1;
%     [a,b] = histc(phot2_id,unique(phot2_id)); y2 = a(b); y2=y2>1;
%     data_photdet1cc(y1)=[]; data_photdet2cc(y2)=[];
%     datadouble12_1cc(y1)=[]; datadouble12_2cc(y2)=[]
% sameTrig=logical(mod(y1+y2,2));
% photDet1(sameTrig,:)=[];
% photDet2(sameTrig,:)=[];

end


%% Get values of accidentals

dataAcc=dataCoinc; dataAcc([data_photdet1,data_photdet2],:)=[];
nAcc=11;ACC_i=zeros(1,nAcc);
plotHist=0;plotNumEvents=0;

% Data structure of allCorrsWithPos
% ind 1=1 for negative delays, ind 1=2 for positive delay, ind 2=1:nAcc is
% number of pulses away. Each cell contains a 2 by (n by 3) matrix
% corresponding to the n counts with the three data information (det,
% clock, time) for det 1 and 2
allCorrDataWithPos=cell(2,nAcc);
hist12_eachAcc=cell(2,nAcc);
allCorrsCount=nan(2,nAcc);

for i=1:nAcc
    pulseTrigDelay_iacc=pulseTrigDelay+i; % Search accidentals at +/- this trig delay
[dataAcc_1,dataAcc_2,dataAcc_photdet1,dataAcc_photdet2]=getDoubleCoincUpdated(pulseTrigDelay_iacc,dataAcc,plotHist,plotNumEvents);

%% Select only counts within correlation window
alltDiffs_iacc=dataAcc_1(:,3)-dataAcc_2(:,3)+2^15;
% Plot below if you want
hisAcc=histcounts(alltDiffs_iacc,0:2^16);plot(timeCoinc,hisAcc); hold on; plot(timeCoinc(corrWind),hisAcc(corrWind))% if you want to plot...
corrCounts_iac=ismember(alltDiffs_iacc,corrWind); %find correlated counts that are within the correlated window corrWind
ACC_i(i)=sum(corrCounts_iac)/2;

dataAcc_1=dataAcc_1(corrCounts_iac,:);dataAcc_2=dataAcc_2(corrCounts_iac,:);
dataAcc_photdet1=dataAcc_photdet1(corrCounts_iac);dataAcc_photdet2=dataAcc_photdet2(corrCounts_iac);

% corrCounts
 %% Can uncomment if you want to keep them... I just comment these for now to save memory
% alldataAcc_1{i}=dataAcc_1;alldataAcc_2{i}=dataAcc_2;
% alldataAcc_photdet1{i}=dataAcc_photdet1;alldataAcc_photdet2{i}=dataAcc_photdet2;
if i==1
    hist12_psAll_iac=histcounts(dataAcc_1(:,3)-dataAcc_2(:,3)+2^15,0:2^16); % keep one just for plotting
end

% Find trig difference to see if pulse corresponds to +/- of pulseTrigDelay_iacc
trigDiffs=dataAcc_1(:,2)-dataAcc_2(:,2); 
trigDiffNeg=trigDiffs<0; trigDiffPos=trigDiffs>0;
allCorrDataWithPos{1,i}=cat(3,dataAcc_1(trigDiffNeg,:),dataAcc_2(trigDiffNeg,:)); allCorrDataWithPos{2,i}=cat(3,dataAcc_1(trigDiffPos,:),dataAcc_2(trigDiffPos,:)); 
allCorrsCount(1,i)=sum(trigDiffNeg); allCorrsCount(2,i)=sum(trigDiffPos); 

% Assing the hist12 difference for plotting the g2 later
hist12_eachAcc{1,nAcc}=histcounts(dataAcc_1(trigDiffNeg,3)-dataAcc_2(trigDiffNeg,3)+2^15,0:2^16); % keep one just for plotting
hist12_eachAcc{2,nAcc}=histcounts(dataAcc_1(trigDiffPos,3)-dataAcc_2(trigDiffPos,3)+2^15,0:2^16); % keep one just for plotting


%%  check if the accidentals should be removed
% % Remove counts found
 dataAcc([dataAcc_photdet1,dataAcc_photdet2],:)=[];
 
 

end
CC=sum(datadouble12_1cc(:,1)); AC=mean(mean(allCorrsCount));mean(ACC_i); 
CAR=(CC-AC)/AC
dCAR=CAR*sqrt(1/CC+1/AC);

figure;bar([fliplr(allCorrsCount(1,:)),CC,allCorrsCount(2,:)]);
title(['g^2_{si} as a function of pulse separation (CAR=' num2str(CAR) '\pm' num2str(dCAR) ')']);
xlabel('Pulse separation'); ylabel('correlated counts')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get all counts and calculate metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n1=sum(hist1_tm) % Per channel, total number of counts after preselection
n2=sum(hist2_tm)


n1Raw=sum(hist1); % Per channel, total number of raw counts
n2Raw=sum(hist2);

c12=numel(datadouble12_1(:,1)) % coincidences between the same beam splitter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate metrics         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% heralding eff and mean photon

% g2_12_34=c12*Np/(n1*n1)

etas=c12/n1
detai=etas*sqrt(1/c12+1/n1)
etai=c12/n2
detas=etai*sqrt(1/c12+1/n2)
% 
% BrightUnNormed=(n1+n2)*(n3+n4)/(c12*totalIntegrationTime) % Pairs/second per mode
% dBrightUnNormed=BrightUnNormed*sqrt(1/(n1+n2)+1/(n3+n4)+1/c12)

BrightUnNormed=(n1)*(n2)/(c12*totalIntegrationTime) % Pairs/second per mode
dBrightUnNormed=BrightUnNormed*sqrt(1/(n1)+1/(n2)+1/c12)
mpn=(n1)*(n2)/(c12*Np)
dmpn=mpn*sqrt(1/(n1)+1/(n2)+1/c12)




h1=figure;
plot(t/1e3,hist1)
hold on
plot(t/1e3,hist2)
% plot(t/1e3,hist3)
% plot(t/1e3,hist4)
plot(t/1e3,hist1_tm)
plot(t/1e3,hist2_tm)
% plot(t/1e3,hist3_tm)
% plot(t/1e3,hist4_tm)
xlabel('Time (ns)');ylabel('Counts')
set(findall(gcf,'-property','FontSize'),'FontSize',16)
%     end
% end
% savefig(h1,[dirDat runDir 'fig-' basefilename vStr pmStr '-' DSP  '.fig']);


h2=figure;
subplot(2,1,1)
plot(timeCoinc,hist12_psAll)
hold on; plot(timeCoinc,hist12_psAll)
% plot(timeCoinc,hist23_psAll); plot(timeCoinc,hist24_psAll)
% xlabel('time (ps)');ylabel('G_2(\tau) ');legend('ch1-ch3','ch1-ch4','ch2-ch3','ch2-ch4')
title(['CAR: ' num2str(CAR)])
subplot(2,1,2)
plot(timeCoinc,hist12_psAll_iac)
hold on; plot(timeCoinc,hist12_psAll_iac)
% plot(timeCoinc,hist23_psAll_iac); plot(timeCoinc,hist24_psAll_iac)
xlabel('time (ps)');ylabel('G_2(\tau) (accidentals)');legend('ch1-ch3','ch1-ch4','ch2-ch3','ch2-ch4')
%         figure;
%         plot(timeCoinc/1e3,hist12_psAll)
%         xlabel('Time (n2)');ylabel('Counts')
%         set(findall(gcf,'-property','FontSize'),'FontSize',16)
% savefig(h2,[dirDat runDir 'fig-' basefilename vStr pmStr '-' DSP  'CoincDelays.fig']);
% 
% close(h1);
% close(h2);








function [coincidenceCount,accidentalsCount,hist13_psAll,hist13_psAll_iac,datadouble13_1,datadouble13_2]=getCarHistCoinc(dataCoinc,nac,...
    pulseTrigDelay,plotHist,plotNumEvents)%postSel,taiInds1,taiInds2,plotHist,plotNumEvents)
%    c1=1;c2=3;
%         dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis

%         coincidenceCount=zeros(1,4); accidentalsCount=zeros(1,4);

%         nac=11;
acc=zeros(1,nac);
for iac=1:nac
    pulseTrigDelay_iac=iac-1+pulseTrigDelay;
    [datadouble13_1_iac,datadouble13_2_iac]=findDoubleCoincSelectedDat(pulseTrigDelay_iac,dataCoinc,plotHist,plotNumEvents);

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
inds1 = ismember(hist1_inds,tm1);                       % Take only the counts there are within the preselection2
data1_tm=data1(inds1,:);
inds1_tm=indsDat1(inds1);
[hist1_tm,~]=histcounts(data1_tm(:,3),0:2^15);
end



function [occurence]=findReptitionNumber(y)
% example:
%y=[ 1 2 3 4 4 4 5 5 6]
% occurence = [ 1     1     1     3     2     1]
f=un1que(y);
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