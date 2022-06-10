function [datadouble_1,datadouble_2,phot2_det1,phot2_det2]=getDoubleCoincUpdated(pulseTrigDelay,dataCoinc,plotHist,plotNumEvents)
%% Find double coincidence
% function is as follows:
% [datadouble12_1,datadouble12_2]=findDoubleCoinc(pulseTrigDelay,1,2,dataSelected,plotHist,plotNumEvents);
% [outputDataCoincidence_fromChannel_c1,outputDataCoincidence_fromChannel_c2]=...
% findDoubleCoinc(Number of pulse trig delays that are recognized as a coincidence,...
% first channel c1, second channel c2 Data under analysis ,plotHist=1 to plot histogram of trigg delays
% ,plotNumEvents=1 to plot the number of coincidences found at each shift);

phot2_det1=[];phot2_det2=[]; 
searchRange=pulseTrigDelay+25; %large search range; stop loop if numEvents goes to 0 after some time.



% dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2),:); % Select only data under analysis

for i=1:searchRange
    
    trigDelays=dataCoinc(:,2)-circshift(dataCoinc(:,2),-i); % Find different between trig delays
    
    if plotHist %% Plot histograms of TrigDelay differences. Useful for debugging.
        figure             
        histogram(trigDelays(abs(trigDelays)<20))
    end
    
    phot_id = [find(trigDelays==-pulseTrigDelay)]; %same pulse trigger
    
    %% needs to be fixed
    numEvents(i)=numel(phot_id); % count number of events found. Use at least 3 shifts; if numEvents fall, exit function.
    if i>2
    if numEvents(i)<sum(numEvents(1:(i-1))) && numEvents(i)==0
        [ 'reached no events after ' num2str(i) ' shifts'];
       break
    end
    end
    
    %% Assign the channel id of each coincidence
    phot_ids_tD2=[phot_id';phot_id'+i]; % Assign the channel id of each coincidence
    
    phot_ids_tD2_Channels=[dataCoinc(phot_ids_tD2(1,:),1)';... % Get Channel
        dataCoinc(phot_ids_tD2(2,:),1)'];
    
    doubleDiffChan= (phot_ids_tD2_Channels(1,:)~=phot_ids_tD2_Channels(2,:)) ; % Only keep coincidences from different channels
    phot_ids_tD2_diffCh= phot_ids_tD2(:,doubleDiffChan);
    phot_ids_tD2_diffCh_Channels=phot_ids_tD2_Channels(:,doubleDiffChan);
    
    [~,phot_ids_tD2_Channels_diffCh_order]=... % Sort counts by channel
        sort(phot_ids_tD2_diffCh_Channels);
    [m,n]=size(phot_ids_tD2_diffCh);
    idx2=repmat(1:n,m,1);
    sortidx=sub2ind([m n],phot_ids_tD2_Channels_diffCh_order,idx2);
    phot_ids_tD2_diffCh_sorted=phot_ids_tD2_diffCh(sortidx);
    
    phot2_det1=[ phot2_det1,phot_ids_tD2_diffCh_sorted(1,:)]; % Assign counts to corresponding channel
    phot2_det2=[ phot2_det2,phot_ids_tD2_diffCh_sorted(2,:)];
    
    
end
    if plotNumEvents %% Useful for debuggin
        figure             
        plot(numEvents)
    end
       datadouble_1=dataCoinc(phot2_det1,:);
       datadouble_2=dataCoinc(phot2_det2,:);

       
       
end
