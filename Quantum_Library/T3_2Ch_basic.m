addpath('/Users/ben/Documents/MATLAB/Instrumentation/Lab Interfacing/HydraHarp_DataTreatment')
addpath( '/Users/ben/Documents/MATLAB/library' )

% Directory with all the data
dirDat='';


runDir='';
basefilename='LargeGap_SigIdlFiltered_68-85_290-310_28dBm_0';
pulseTrigDelay=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Script %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




load([basefilename '.mat']);% dirDat basefilename  '_' num2str(fileNums-1) '.mat'])
data(data(:,1)==0,:)=[];
lastGarbageData=find(data(1:3e3,1)>4,1,'last');
data(1:(lastGarbageData+1),:)=[];

%%% This section can be used to reassign certain channels
% % data(data(:,1)==1,:)=[];
%  data(data(:,1)==2,1)=3;
% data(data(:,1)==3,1)=2;
% data(data(:,1)==4,1)=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic Data Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This section does basic processing namely assign1ng data per
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

res=256; % HHresolution. Ok to hard code (ps)
t = (1:2^15)*res;
lastInd=find(hist1>5,1,'last'); % Last index of the histogram with mean1nfull information
trigPeriod=t(lastInd)*1e-12; % Correponds to period of pulsed signal
totalIntegrationTime=(data(end,2)-data(1,2))*trigPeriod;
Np=data(end,2)-data(1,2); % Number of pulses is simply the total number of trigger clocks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-selection Data Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

histPks=zeros(4,1);


preSelectionWidth=15; 
range=-preSelectionWidth:preSelectionWidth;
[maxVal1,maxInd1]=max(hist1); time_win{1}=maxInd1+range;
[maxVal2,maxInd2]=max(hist2); time_win{2}=maxInd2+range;

[inds1,data1_tm,inds1_tm,hist1_tm]=preSelect(time_win{1},hist1_inds,data1,indsDat1);
[inds2,data2_tm,inds2_tm,hist2_tm]=preSelect(time_win{2},hist2_inds,data2,indsDat2);

[~,order]=sort([inds1_tm,inds2_tm]);%,inds3_tm,inds4_tm]);
dataSelected=[data1_tm;data2_tm];%;data3_tm;data4_tm]; 
dataSelected=dataSelected(order,:);   % Put it back into the order at which it was detected for the "trigger selection"

%% Plot selected region2
%
figure
plot(hist1);hold on; plot(hist1_tm)
plot(hist2);hold on; plot(hist2_tm)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Coincidences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find double and triple coincidences %%%%%%%%%%%%%%%%%%%%%%%%%
% function is as follows:
% [datadouble12_1,datadouble12_2]=findDoubleCoinc(pulseTrigDelay,1,2,dataSelected,plotHist,plotNumEvents);
% maps to:
% [outputDataCoincidence_fromChannel_c1,outputDataCoincidence_fromChannel_c2]=...
% findDoubleCoinc(Number of pulse trig delays that are recogn1zed as a coincidence,...
% first channel label for c1, second channel c2 ,plotHist=1 to plot histogram of trigg delays
% ,plotNumEvents=1 to plot the number of coincidences found at each
% shift); Similarly for triple coincidence function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



plotHist=1;plotNumEvents=1;
nac=11;

c1=1;c2=2;

dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis

[datadouble12_1,datadouble12_2]=findDoubleCoincSelectedDat(pulseTrigDelay,dataCoinc,plotHist,plotNumEvents);

plotHist=0;plotNumEvents=0;
[coincidenceCount,accidentalsCount,hist14_psAll,hist14_psAll_iac,datadouble12_1,datadouble12_2]=getCarHistCoinc(dataCoinc,nac,...
    pulseTrigDelay,plotHist,plotNumEvents);
CC=sum(coincidenceCount); AC=sum(accidentalsCount); CAR=(CC-AC)/AC; dCAR=CAR*sqrt(1/CC+1/AC);

[hist12_psAll,~,hist12_psAll_inds]=histcounts(datadouble12_1(:,3)-datadouble12_2(:,3)+2^15,0:2^16);
timeCoinc = ((1:2^16)-2^15)*res;

figure;plot(timeCoinc,hist12_psAll);


%% Plot correlated counts on t1-t2 axis





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get all counts and calculate metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

BrightUnNormed=(n1+n2)*(n3+n4)/(c12*totalIntegrationTime) % Pairs/second per mode
dBrightUnNormed=BrightUnNormed*sqrt(1/(n1+n2)+1/(n3+n4)+1/c12)
mpn=(n1)*(n2)/(c12*Np)
dmpn=mpn*sqrt(1/(n1+n2)+1/(n3+n4)+1/c12)


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
    'c12','n1','n2','Np');
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
xlabel('Time (n2)');ylabel('Counts')
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
%         xlabel('Time (n2)');ylabel('Counts')
%         set(findall(gcf,'-property','FontSize'),'FontSize',16)
savefig(h2,[dirDat runDir 'fig-' basefilename vStr pmStr '-' DSP  'CoincDelays.fig']);

close(h1);
close(h2);





%
%     if iTAI==1
%         newFig=1;hetai=[]; hetas=[];hmpn=[]; hBright=[];
%         hnRaw=[]; hnCh=[]; hcoinc=[]; hCAR=[];
%         saveIt=0;  hg2h_i_vi=[]; hg2h_s_vi=[];
%         hg2_12_vi=[]; hg2_34_vi=[];
%     elseif iTAI==numTaiDSP
%         newFig=0; saveIt=1; % On last pass, save figure;
%     else
%         newFig=0; saveIt=0;
%     end
%
%     hCAR=plotMetric(pows,CAR_vi,dCAR_vi,newFig,hCAR,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%     hetas=plotMetric(pows,etas_vi,detas_vi,newFig,hetas,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%     hetai=plotMetric(pows,etai_vi,detai_vi,newFig,hetai,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%     hmpn=plotMetric(pows,mpn_vi,dmpn_vi,newFig,hmpn,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%     NormalizedBrightness=BrightUnNormed_vi./(pows*tapLossFac);
%     dNormalizedBrightness=dBrightUnNormed_vi./(pows*tapLossFac);
%     hBright=plotMetric(pows,NormalizedBrightness,dNormalizedBrightness,newFig,hBright,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%
%     hnRaw=plotMetric(pows,nRaw_vi,sqrt(nRaw_vi),newFig,hnRaw,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%     hnCh=plotMetric(pows,nCh_vi,(sqrt(nCh_vi)),newFig,hnCh,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%     hcoinc=plotMetric(pows,coinc_vi,sqrt(coinc_vi),newFig,hcoinc,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%
%     hg2h_i_vi=plotMetric(pows,g2h_i_vi,dg2h_i_vi,newFig,hg2h_i_vi,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%     hg2h_s_vi=plotMetric(pows,g2h_s_vi,dg2h_s_vi,newFig,hg2h_s_vi,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%     hg2_12_vi=plotMetric(pows,g2_12_vi,dg2_12_vi,newFig,hg2_12_vi,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%     hg2_34_vi=plotMetric(pows,g2_34_vi,dg2_34_vi,newFig,hg2_34_vi,saveIt,[pmStr, '-' DSP],[dirDat runDir],basefilename);
%
%
%     nCh_vi(vi,:)=[n1,n2,n3,n4];
%     coinc_vi(vi,:)=[c13,c14,c23,c24];
%     selfCoinc_vi(vi,:)=[c12,c34];
%     triples_vi(vi,:)=[c123,c124,c134,c234];
%     nRaw_vi(vi,:)=[n1Raw,n2Raw,n3Raw,n4Raw];
%








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