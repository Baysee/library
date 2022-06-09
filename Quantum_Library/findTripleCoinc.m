function [dataTrip_1,dataTrip_2,dataTrip_3]=findTripleCoinc(pulseTrigDelay,c1,c2,c3,dataSelected,plotNumEvents)
%% Find triple coincidence
% Right now, this code does a first search by finding differences between
% each element and the element 1 to "search range" away. If the difference
% returns  a trig delay difference equal to "pulseTrighDelay", this means that they are
% coincident. For each coincidence, look if next element also occured
% during same trig delay (i.e., check if next element is also on the same trigger, assuming the pulse Trig Delay is 0)
%. If so, then identify as triple coincidence, and
% keep if the channels are all different. Next change would be to do a
% second scanning range on the third event; this will be necessary if the
% trig delay identified for coincidence is not 0.

phot3_det1=[];phot3_det2=[]; phot3_det3=[];photDet4=[];

searchRange=pulseTrigDelay+5;

dataCoinc=dataSelected( (dataSelected(:,1)==c1) | (dataSelected(:,1)==c2)| (dataSelected(:,1)==c3),:);

nCh=3;
for i=1:searchRange
    
   trigDelays=dataCoinc(:,2)-circshift(dataCoinc(:,2),-i); % Find different between trig delays
   tD=(trigDelays==pulseTrigDelay); % find differences equal to definition of pulseTrigDelay corresponsding to a coincident (usually 0)
   tD3=tD+circshift(tD,-1); % Sum logical array to find two subsequent successfull coincidence (means a triple)
   phot_id_tD3=find(tD3==2);
   
   % check if events are still being found
   numEvents(i)=numel(phot_id_tD3);
    if i>1
    if numEvents(i)<numEvents(i-1) && numEvents(i)==0
       [ 'reached no events after ' num2str(i) ' shifts']
       break
    end
    end
   
   phot_ids_tD3=[phot_id_tD3';phot_id_tD3'+1;phot_id_tD3'+2]; % Assign the channel id of each triple
   phot_ids_tD3_Channels=[dataCoinc(phot_ids_tD3(1,:),1)';...
       dataCoinc(phot_ids_tD3(2,:),1)';dataCoinc(phot_ids_tD3(3,:),1)'];

     tripleDiffChan= (phot_ids_tD3_Channels(1,:)~=phot_ids_tD3_Channels(2,:)) & ... % take counts from different channels only
       (phot_ids_tD3_Channels(1,:)~=phot_ids_tD3_Channels(3,:)) & (phot_ids_tD3_Channels(2,:)~=phot_ids_tD3_Channels(3,:));
    phot_ids_tD3_diffCh= phot_ids_tD3(:,tripleDiffChan);
    phot_ids_tD3_diffCh_Channels=phot_ids_tD3_Channels(:,tripleDiffChan);
     [phot_ids_tD3_diffCh_Channels_sorted,phot_ids_tD3_Channels_diffCh_order]=sort(phot_ids_tD3_diffCh_Channels);
%    phot_ids_tD3=phot_ids_tD3';
    [m,n]=size(phot_ids_tD3_diffCh);
    idx2=repmat(1:n,m,1);
    sortidx=sub2ind([m n],phot_ids_tD3_Channels_diffCh_order,idx2);
       phot_ids_tD3_diffCh_sorted=phot_ids_tD3_diffCh(sortidx);

       phot3_det1=[ phot3_det1,phot_ids_tD3_diffCh_sorted(1,:)];
       phot3_det2=[ phot3_det2,phot_ids_tD3_diffCh_sorted(2,:)];
       phot3_det3=[ phot3_det3,phot_ids_tD3_diffCh_sorted(3,:)];
   
end
    if plotNumEvents %% Useful for debuggin
        figure             
        plot(numEvents)
    end
       dataTrip_1=dataCoinc(phot3_det1,:);
       dataTrip_2=dataCoinc(phot3_det2,:);
       dataTrip_3=dataCoinc(phot3_det3,:);
end